#!/usr/bin/env python
"""Build a TSV of read_id, query_transcript, ref_transcript gene."""
import argparse

import pandas as pd


def parse_args():
    """Create argument parser."""
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "--read_tr_map",
        help="TSV of read to reference transcript name mapping"
    )

    parser.add_argument(
        "--gffcompare_tmap",
        help="The .tmap output from gffcompare"
    )

    parser.add_argument(
        "--all_read_ids",
        help="Ordered read ids of mapped and unampped reads"
    )

    parser.add_argument(
        "--output",
        help="Path to output file"
    )
    return parser.parse_args()


def main(read_tr_map, tmap, read_order, outfile):
    """Entry point.

    read_tr_map: maps reads to query transcripts
    tamp: maps query isoforms to reference isoforms
    """
    df_r = pd.read_csv(
        read_tr_map,
        names=['read_id', 'qry_id'],
        index_col=None,
        sep='\t')
    df_t = pd.read_csv(
        tmap,
        index_col=None,
        sep='\t')

    # Merge gffcompare annotation to read_id
    df_tr = df_r.merge(
        df_t[['qry_id', 'ref_id', 'ref_gene_id', 'class_code']],
        left_on='qry_id',
        right_on='qry_id')

    df_read_order = pd.read_csv(read_order, index_col=None, names=['read_id'])
    df_tr['status'] = 'Assigned'
    df_out = df_tr.merge(
        df_read_order, left_on='read_id', right_on='read_id', how='right')
    df_out.loc[df_out['status'] != 'Assigned', 'status'] = 'Unasssigned'

    # sort out read duplicated
    df_out['dup'] = df_out['read_id'].duplicated(keep=False)
    df_dups = df_out[df_out['dup']]
    df_out = df_out.loc[~df_out['dup']]

    # Duplicate reads occur when Salmon maps anbiquously to multiple query
    # transcritps. If all map to same reference transcript, assign to it,
    # else assign as ambiguous.
    # TODO: Can we get a mapping score from Salmon to help here.
    # Should a class of '=' automatically trump other codes?
    # Also are there other reasons we get duplicates?
    results = []
    for read_id, dup_frame in df_dups.groupby('read_id'):
        if all(dup_frame['ref_id'].values[0] == dup_frame['ref_id']):
            results.append(dup_frame.iloc[0])
        else:
            row = dup_frame.iloc[1]
            row['status'] = 'Unassigned ambiguous'
            results.append(row)

    resolved_dups = pd.DataFrame.from_records(results)
    df_out = pd.concat([df_out, resolved_dups])
    df_out.set_index('read_id', drop=True, inplace=True)

    # set status based on gffcompare class code
    df_out.loc[df_out['class_code'].isin(['i', 'p', 's', 'u']), 'status'] = \
        'Unassigned gffcompare class'

    df_out.drop(columns=['dup'], inplace=True)
    df_out = df_out.loc[df_read_order['read_id']]
    df_out.fillna('-', inplace=True)
    df_out.to_csv(outfile, sep='\t')


if __name__ == "__main__":
    args = parse_args()

    main(
        args.read_tr_map, args.gffcompare_tmap, args.all_read_ids, args.output)
