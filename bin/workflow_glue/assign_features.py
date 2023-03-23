#!/usr/bin/env python
"""Assign reference genes and transcripts to read_id."""

import re

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("assign_features")

    parser.add_argument(
        "--query_transcript_read_assign",
        help="TSV with columns: read_id, query_transcript"
    )

    parser.add_argument(
        "--gffcompare_tmap",
        help="The .tmap output from gffcompare"
    )

    parser.add_argument(
        "--gtf",
        help="GTF for chromosome. Should contain gene_id and gene_name "
    )

    parser.add_argument(
        "--tags",
        help="TSV file with read_id index and genomic alignment mapq column."
    )

    parser.add_argument(
        "--output",
        help="Path to output file"
    )

    parser.add_argument(
        "--min_mapq",
        help="mapq threshold from genomic alignment above which reads will be assigned "
             "a gene and transcript",
        type=int, default=30
    )
    return parser


def parse_gtf(annotation_file):
    """Get a mapping of query transcript id to reference gene name."""
    with open(annotation_file, 'r') as fh:

        records = []
        for line in fh:
            if line.split('\t')[2] != 'transcript':
                continue
            g = re.search(r'gene_name "(.+?)"', line)
            if g:
                gene_name = g.group(1)
            else:
                gene_name = None
            t = re.search(r'transcript_id "(.+?)"', line)
            if t:
                transcript_id = t.group(1)
            else:
                transcript_id = None

            records.append([transcript_id, gene_name])
        df = pd.DataFrame.from_records(records, columns=['ref_id', 'gene_name'])
        return df


def main(args):
    """Entry point.

    Merge the stringtie output containing read_id and query_transcript (transcript built
    by Stringtie) with the gffcompare output, which maps query transcript to reference
    transcritpt.
    """
    logger = get_named_logger('AssignFeat')
    logger.info('Mapping reference info to reads.')

    # Mapq is from the genome alignemnt
    df_mapq = pd.read_csv(args.tags, sep='\t', index_col=0, usecols=['read_id', 'mapq'])

    # Dataframe with read_id and query_transcript (transcript built by strintie)
    df_query_transcript = pd.read_csv(
        args.query_transcript_read_assign,
        index_col=None,
        sep='\t')

    # DataFrame which maps query transcript to reference transcript
    df_gffcompare_tmap = pd.read_csv(
        args.gffcompare_tmap,
        index_col=None,
        sep='\t')

    # Merge on query_id to assign read_id to reference transcript.
    # Only keep reads that have been assigned a query_id
    df_tr = df_query_transcript.merge(
        df_gffcompare_tmap[['qry_id', 'ref_id', 'class_code']],
        left_on='qry_id',
        right_on='qry_id')

    # Merge the gene names from the original input annotation gtf.
    # It is possible to get the gene name and reference transcript map from the
    # gffcompare *.refmap file, but this file seems to be been missing
    # some query_ids/ref_ids.
    df_ann = parse_gtf(args.gtf)
    df_tr = df_tr.merge(df_ann, how='left', left_on='ref_id', right_on='ref_id')
    df_tr.gene_name = df_tr.gene_name.fillna('-')
    df_tr.rename(
        columns={
            'ref_id': 'transcript',
            'gene_name': 'gene',
        }, inplace=True)

    # Set status based on gffcompare class code
    # If any of these classes, call transcript as unclassified.
    df_tr.loc[df_tr['class_code'].isin(['i', 'p', 's', 'u']), 'transcript'] = '-'

    # If genomic mapq below threshold, set both gene and transcript status to unassigned
    df_tr = df_tr.merge(df_mapq, how='left', left_on='read_id', right_on='read_id')
    df_tr.loc[df_tr.mapq < args.min_mapq, 'gene'] = '-'
    df_tr.loc[df_tr.mapq < args.min_mapq, 'transcript'] = '-'

    df_tr.set_index('read_id', drop=True, inplace=True)

    df_tr.fillna('-', inplace=True)
    df_tr.to_csv(args.output, sep='\t')
