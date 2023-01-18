#!/usr/bin/env python
"""Cell Umi gene stats."""
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("UmiGene")

    # Positional mandatory arguments
    parser.add_argument(
        "tsv",
        help="TSV file with read ID, cell, UMI, and gene",
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output TSV file containing table of read ID, cell barcode, UMI, \
        and gene values extracted from the input BAM [cell_umi_gene.tsv]",
        default="cell_umi_gene.tsv",
    )

    return parser


def read_table(args):
    """Read table."""
    df = pd.read_csv(args.tsv, sep="\t")
    return df


def per_cell_stats(df):
    """Run groupby on cell barcode and return table of per-cell stats."""
    reads_per_cell = df.groupby("barcode")["read_id"].count()
    umis_per_cell = df.groupby("barcode")["umi"].nunique()
    genes_per_cell = df.groupby("barcode")["gene"].nunique()

    return pd.concat([reads_per_cell, umis_per_cell, genes_per_cell], axis=1)


def per_cell_umi_stats(df):
    """Run groupby on cell+UMI and return table of per-cell+UMI stats."""
    reads_per_cell_umi = df.groupby(["barcode", "umi"])["read_id"].count()
    genes_per_cell_umi = df.groupby(["barcode", "umi"])["gene"].nunique()

    return pd.concat([reads_per_cell_umi, genes_per_cell_umi], axis=1)


def per_cell_gene_stats(df):
    """Run groupby on cell+gene and return table of per-cell+gene stats."""
    reads_per_cell_gene = df.groupby(["barcode", "gene"])["read_id"].count()
    umis_per_cell_gene = df.groupby(["barcode", "gene"])["umi"].nunique()

    return pd.concat([reads_per_cell_gene, umis_per_cell_gene], axis=1)


def per_gene_bulk_stats(df):
    """Run groupby on gene and return table of per-gene (bulk) stats."""
    reads_per_gene_bulk = df.groupby("gene")["read_id"].count()
    umis_per_gene_bulk = df.groupby("gene")["umi"].nunique()

    return pd.concat([reads_per_gene_bulk, umis_per_gene_bulk], axis=1)


def main(args):
    """Run entry point."""
    logger = get_named_logger('cell_umi_gene_stats')

    df = read_table(args)

    n_igk_reads = df.loc[df["gene"].str.startswith("IGK"), :].shape[0]
    n_igl_reads = df.loc[df["gene"].str.startswith("IGL"), :].shape[0]
    n_igh_reads = df.loc[df["gene"].str.startswith("IGH"), :].shape[0]

    logger.info(f"IGH reads: {n_igh_reads}")
    logger.info(f"IGK reads: {n_igk_reads}")
    logger.info(f"IGL reads: {n_igl_reads}")


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
