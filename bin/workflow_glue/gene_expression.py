#!/usr/bin/env python
"""Gene expression."""
from pathlib import Path

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("gene_expression")

    # Positional mandatory arguments
    parser.add_argument(
        "--read_tags",
        help="TSV with read_id and associated tags.",
        type=Path)

    # Optional arguments
    parser.add_argument(
        "--output_prefix",
        help="Output TSV file containing gene expression matrix, where genes \
        are rows and barcodes are columns [gene_expression.tsv]",
        default="gene_expression.tsv",
    )

    return parser


def process_tag_tsv(read_tags_tsv):
    """Convert TSV of read data to gene and transcript expression matrices.

    :param read_tags_tsv: read_tags_tsv.
    :type read_tags_tsv: Path
    """
    df = pd.read_csv(read_tags_tsv, sep='\t', index_col=0)

    def process_dataframe(feature='gene'):
        dfg = df[[feature, 'CB', 'UB']]
        dfg = dfg.groupby([feature, 'CB']).nunique()['UB'].reset_index()
        dfg = pd.pivot(dfg, index=feature, columns='CB', values='UB')
        dfg = dfg.dropna(how='all', axis=1)
        dfg.fillna(0, inplace=True)
        dfg.drop('-', axis=1, inplace=True, errors='ignore')  # Drop unassigned cells
        dfg.drop('-', axis=0, inplace=True, errors='ignore')  # Drop unassigned genes
        return dfg

    df_gene = process_dataframe('gene')
    df_transcript = process_dataframe('transcript')
    df_transcript = df_transcript.drop('-', errors='ignore')

    return df_gene, df_transcript


def main(args):
    """
    Iterate through the BAM file and count unique UMIs (UB tag) \
    associated with each gene (GN tag) and cell barcode (CB tag).

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    logger = get_named_logger('GeneEx')
    logger.info(
        f"Building gene/transcript expression matrices from {args.read_tags}")

    gene_expression_df, transcript_expression_df = \
        process_tag_tsv(args.read_tags)

    gene_expression_df.to_csv(
        f'{args.output_prefix}.gene_expression.counts.tsv',
        sep="\t", index_label="gene")
    transcript_expression_df.to_csv(
        f'{args.output_prefix}.transcript_expression.counts.tsv',
        sep="\t", index_label="transcript")


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
