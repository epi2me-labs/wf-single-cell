import argparse
import logging
import os
import sys

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "counts",
        help="Matrix of read counts per gene (row) \
        per cell (column)",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--min_genes",
        help="Filter out cells that contain fewer than <min_genes> genes [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--min_cells",
        help="Filter out genes that are observed in fewer than <min_cells> \
        cells [3]",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--max_mito",
        help="Filter out cells where more than <max_mito> percent of counts \
        belong to mitochondrial genes [5]",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--norm_count",
        help="Normalize to this number of counts per cell [10000]",
        type=int,
        default=10000,
    )

    parser.add_argument(
        "--output",
        help="Output TSV file containing processed gene expression matrix, \
            where genes are rows and cells are columns [expression.processed.tsv]",
        type=str,
        default="expression.processed.tsv",
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


def filter_cells(df, args):
    """
    Remove cells that express fewer than N=<args.min_genes> unique genes.
    """
    df = df.transpose()
    df["total"] = df.sum(axis=1)

    # Remove cells that have fewer than <args.min_genes> unique genes
    df["n_genes"] = df.astype(bool).sum(axis=1)
    n_dropped = df[df["n_genes"] < args.min_genes].shape[0]
    logger.info(f"Dropping {n_dropped} cells with < {args.min_genes} genes")
    df = df[df["n_genes"] >= args.min_genes]
    df = df.drop(["n_genes"], axis=1)

    # Filter out cells where mitochondrial genes comprise more than
    # <args.max_mito> percentage of the total count
    mito_genes = [gene for gene in df.columns if gene.find("MT-") == 0]
    df["mito_total"] = df.loc[:, mito_genes].sum(axis=1)
    df["mito_pct"] = 100 * df["mito_total"] / df["total"]
    n_mito = df[df["mito_pct"] > args.min_genes].shape[0]
    logger.info(f"Dropping {n_mito} cells with > {args.max_mito}% mitochondrial reads")
    df = df[df["mito_pct"] <= args.max_mito]

    df = df.drop(["mito_total", "mito_pct"], axis=1)
    df = df.transpose()

    if df.shape[1] == 0:
        raise Exception("All cells have been filtered out!")

    return df


def filter_genes(df, args):
    """ """
    df["n_cells"] = df.astype(bool).sum(axis=1)

    # Remove genes that are present in fewer than <args.min_cells> unique cells
    n_dropped = df[df["n_cells"] < args.min_cells].shape[0]
    logger.info(f"Dropping {n_dropped} genes observed in < {args.min_cells} cells")
    df = df[df["n_cells"] >= args.min_cells]
    df = df.drop(["n_cells"], axis=1)

    if df.shape[0] == 0:
        raise Exception("All genes have been filtered out!")

    return df


def normalize(df, args):
    """ """
    df = df.transpose()

    # cell_count / cell_total = X / <args.norm_count>
    logger.info(f"Normalizing counts to {args.norm_count} reads per cell")
    genes = df.columns != "total"
    df.loc[:, genes] = df.loc[:, genes] * args.norm_count
    df.loc[:, genes] = df.loc[:, genes].div(df["total"], axis=0)
    df = df.drop(["total"], axis=1)
    df = df.transpose()

    return df


def logarithmize(df):
    """ """
    # Add pseudo count and logarithmize the normalized counts
    logger.info("Logarithmizing counts")
    df = np.log10(df + 1)
    return df


def main(args):
    init_logger(args)

    df = pd.read_csv(args.counts, sep="\t").set_index("gene")

    df = filter_cells(df, args)
    df = filter_genes(df, args)
    df = normalize(df, args)
    df = logarithmize(df)

    logger.info(f"Processed matrix: {df.shape[0]} genes x {df.shape[1]} cells")
    df.to_csv(args.output, sep="\t")


if __name__ == "__main__":
    args = parse_args()

    main(args)
