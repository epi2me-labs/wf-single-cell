import argparse
import logging
import os
import sys

import pandas as pd
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "tsv",
        help="TSV file with read ID, cell, UMI, and gene",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output TSV file containing table of read ID, cell barcode, UMI, \
        and gene values extracted from the input BAM [cell_umi_gene.tsv]",
        type=str,
        default="cell_umi_gene.tsv",
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
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def read_table(args):
    """ """
    df = pd.read_csv(args.tsv, sep="\t")
    return df


def per_cell_stats(df):
    """
    Run groupby on cell barcode and return table of per-cell stats.
    """
    reads_per_cell = df.groupby("barcode")["read_id"].count()
    umis_per_cell = df.groupby("barcode")["umi"].nunique()
    genes_per_cell = df.groupby("barcode")["gene"].nunique()

    return pd.concat([reads_per_cell, umis_per_cell, genes_per_cell], axis=1)


def per_cell_umi_stats(df):
    """
    Run groupby on cell+UMI and return table of per-cell+UMI stats.
    """
    reads_per_cell_umi = df.groupby(["barcode", "umi"])["read_id"].count()
    genes_per_cell_umi = df.groupby(["barcode", "umi"])["gene"].nunique()

    return pd.concat([reads_per_cell_umi, genes_per_cell_umi], axis=1)


def per_cell_gene_stats(df):
    """
    Run groupby on cell+gene and return table of per-cell+gene stats.
    """
    reads_per_cell_gene = df.groupby(["barcode", "gene"])["read_id"].count()
    umis_per_cell_gene = df.groupby(["barcode", "gene"])["umi"].nunique()

    return pd.concat([reads_per_cell_gene, umis_per_cell_gene], axis=1)


def per_gene_bulk_stats(df):
    """
    Run groupby on gene and return table of per-gene (bulk) stats.
    """
    reads_per_gene_bulk = df.groupby("gene")["read_id"].count()
    umis_per_gene_bulk = df.groupby("gene")["umi"].nunique()

    return pd.concat([reads_per_gene_bulk, umis_per_gene_bulk], axis=1)


def main(args):
    init_logger(args)

    df = read_table(args)

    n_IGK_reads = df.loc[df["gene"].str.startswith("IGK"), :].shape[0]
    n_IGL_reads = df.loc[df["gene"].str.startswith("IGL"), :].shape[0]
    n_IGH_reads = df.loc[df["gene"].str.startswith("IGH"), :].shape[0]

    logger.info(f"IGH reads: {n_IGH_reads}")
    logger.info(f"IGK reads: {n_IGK_reads}")
    logger.info(f"IGL reads: {n_IGL_reads}")

    # df_gene_bulk = per_gene_bulk_stats(df)
    # print(df_gene_bulk.head())
    # print(df_gene_bulk.describe())
    #
    # df_cell = per_cell_stats(df)
    # print(df_cell.head())
    #
    # df_cell_umi = per_cell_umi_stats(df)
    # print(df_cell_umi.head())
    #
    # df_cell_gene = per_cell_gene_stats(df)
    # print(df_cell_gene.head())
    # print(df_cell_gene.describe())


if __name__ == "__main__":
    args = parse_args()

    main(args)
