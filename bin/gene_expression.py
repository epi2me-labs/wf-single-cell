#!/usr/bin/env python
"""Gene expression."""
import argparse
import logging
from pathlib import Path

import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    """Create argument parser."""
    parser = argparse.ArgumentParser()

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

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    """Initiate logger."""
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def process_tag_tsv(read_tags_tsv):
    """Convert TSV of read data to gene and transcript expression matrices.

    :param read_tags_tsv: read_tags_tsv.
    :type read_tags_tsv: Path
    """
    # Build regular expression for "gene" annotations where no gene was found
    regex = r"[a-zA-Z0-9]+_\d+_\d+"  # e.g. chr7_44468000_44469000

    df = pd.read_csv(read_tags_tsv, sep='\t', index_col=0)

    def process_dataframe(feature='gene'):
        dfg = df[[feature, 'barcode', 'umi']]
        dfg = dfg.groupby([feature, 'barcode']).nunique()['umi'].reset_index()
        dfg = pd.pivot(dfg, index=feature, columns='barcode', values='umi')
        dfg = dfg.dropna(how='all', axis=1)
        dfg.fillna(0, inplace=True)
        dfg = dfg.loc[~dfg.index.str.contains(regex, regex=True)]
        return dfg

    df_gene = process_dataframe('gene')
    df_transcript = process_dataframe('transcript')
    df_transcript = df_transcript.drop('-')

    return df_gene, df_transcript


def process_reads(args):
    """
    Iterate through the BAM file and count unique UMIs (UB tag) \
    associated with each gene (GN tag) and cell barcode (CB tag).

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
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


def main(args):
    """Run entry point."""
    init_logger(args)

    process_reads(args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
