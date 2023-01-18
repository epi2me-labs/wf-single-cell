#!/usr/bin/env python
"""Build a TSV of read_id, query_transcript, ref_transcript gene."""

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("isoform_read_mapping")

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
    return parser


def main(args):
    """Entry point.

    read_tr_map: maps reads to query transcripts
    tamp: maps query isoforms to reference isoforms
    """
    logger = get_named_logger('IsoReadMap')
    logger.info('Mapping reads to transcript isoforms')
    df_r = pd.read_csv(
        args.read_tr_map,
        names=['read_id', 'qry_id'],
        index_col=None,
        sep='\t')
    df_t = pd.read_csv(
        args.gffcompare_tmap,
        index_col=None,
        sep='\t')

    # Merge gffcompare annotation to read_id
    df_tr = df_r.merge(
        df_t[['qry_id', 'ref_id', 'ref_gene_id', 'class_code']],
        left_on='qry_id',
        right_on='qry_id')

    df_read_order = pd.read_csv(
        args.all_read_ids, index_col=None, names=['read_id'])
    df_tr['status'] = 'Assigned'
    df_out = df_tr.merge(
        df_read_order, left_on='read_id', right_on='read_id', how='right')
    df_out.loc[df_out['status'] != 'Assigned', 'status'] = 'Unasssigned'

    # sort out read duplicated
    df_out['dup'] = df_out['read_id'].duplicated(keep=False)
    df_dups = df_out[df_out['dup']]
    df_out = df_out.loc[~df_out['dup']]

    # If all map to same reference transcript, assign to it,
    # else assign as ambiguous.
    # Should a class of '=' automatically trump other codes?
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
    df_out.to_csv(args.output, sep='\t')


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
