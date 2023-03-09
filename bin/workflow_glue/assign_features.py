#!/usr/bin/env python
"""Assign reference genes and transcripts to read_id."""

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


def main(args):
    """Entry point.

    Merge the stringtie output containing read_id and query_transcript (transcript built
    by Stringtie) with the gffcompare output, which maps query transcript to reference
    transcritpt.
    """
    logger = get_named_logger('AssignFeat')
    logger.info('Mapping reference info to reads.')

    df_mapq = pd.read_csv(args.tags, sep='\t', index_col=0, usecols=['read_id', 'mapq'])

    # Dataframe with read_id and wuery_transcript (transcript built by strintie)
    df_stringtie_assign = pd.read_csv(  # Add headers to the files
        args.query_transcript_read_assign,
        names=['read_id', 'qry_id'],
        index_col=None,
        sep='\t')

    # DataFrame whcih maps query transcript to reference transcript
    df_gffcompare_tmap = pd.read_csv(
        args.gffcompare_tmap,
        index_col=None,
        sep='\t')

    # Merge on query_id to assign read_id to reference gene and transcript.
    df_tr = df_stringtie_assign.merge(
        df_gffcompare_tmap[['qry_id', 'ref_id', 'ref_gene_id', 'class_code']],
        left_on='qry_id',
        right_on='qry_id')

    df_tr.rename(columns={
        'ref_id': 'transcript',
        'ref_gene_id': 'gene',
    }, inplace=True)

    # set status based on gffcompare class code
    # If any of these classes, call transcript as unclassified
    df_tr.loc[df_tr['class_code'].isin(['i', 'p', 's', 'u']), 'transcript'] = '-'

    # If genomic mapq below threshold, set both gene and transcript status to unassigned
    df_tr = df_tr.merge(df_mapq, left_on='read_id', right_on='read_id')
    df_tr.loc[df_tr.mapq < args.min_mapq, 'gene'] = '-'
    df_tr.loc[df_tr.mapq < args.min_mapq, 'transcript'] = '-'

    df_tr.set_index('read_id', drop=True, inplace=True)

    df_tr.fillna('-', inplace=True)
    df_tr.to_csv(args.output, sep='\t')
