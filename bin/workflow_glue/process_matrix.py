#!/usr/bin/env python
"""Process matrix."""
import numpy as np
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("process_matrix")

    parser.add_argument(
        "--gene_counts",
        help="Matrix of read counts per gene (row) \
        per cell (column)"
    )

    parser.add_argument(
        "--transcript_counts",
        help="Matrix of read counts per transcript (row) \
        per cell (column)"
    )

    # Optional arguments
    parser.add_argument(
        "--min_genes",
        help="Filter out cells that contain fewer than \
        <min_genes> genes [100]",
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
        "--mito_prefix",
        help="prefix(s) to identify mitochondrial genes. Multiple can be \
            supplied as comma-sperated prefixes",
        default="MT-",
    )

    parser.add_argument(
        "--sample_id",
        help="prefix for output file name"
    )

    parser.add_argument(
        "--norm_count",
        help="Normalize to this number of counts per cell [10000]",
        type=int,
        default=10000,
    )

    return parser


def filter_cells(df_gene, df_tr, args):
    """Remove cells that express fewer than N=<args.min_genes> \
        unique genes."""
    logger = get_named_logger('ProcMatrix')
    df_gene = df_gene.transpose()
    df_tr = df_tr.transpose()
    df_gene["total"] = df_gene.sum(axis=1)
    df_tr["total"] = df_tr.sum(axis=1)

    # Remove cells that have fewer than <args.min_genes> unique genes
    df_gene["n_genes"] = df_gene.astype(bool).sum(axis=1)
    n_dropped = df_gene[df_gene["n_genes"] < args.min_genes].shape[0]
    logger.info(f"Dropping {n_dropped} cells with < {args.min_genes} genes")
    df_gene = df_gene[df_gene["n_genes"] >= args.min_genes]
    df_gene = df_gene.drop(["n_genes"], axis=1)

    # Filter out cells where mitochondrial genes comprise more than
    # <args.max_mito> percentage of the total count
    mito_prefixes = args.mito_prefix.strip().split(',')
    mito_genes = []
    for gene in df_gene.columns:
        for mpre in mito_prefixes:
            if gene.find(mpre) == 0:
                mito_genes.append(gene)

    df_gene["mito_total"] = df_gene.loc[:, mito_genes].sum(axis=1)
    df_gene["mito_pct"] = 100 * df_gene["mito_total"] / df_gene["total"]
    df_gene["mito_pct"].to_csv(f"{args.sample_id}.gene_expression.mito.tsv", sep="\t")
    n_mito = df_gene[df_gene["mito_pct"] > args.max_mito].shape[0]
    logger.info(
        f"Dropping {n_mito} cells with > {args.max_mito}% mitochondrial reads")
    df_gene = df_gene[df_gene["mito_pct"] <= args.max_mito]

    df_gene = df_gene.drop(["mito_total", "mito_pct"], axis=1)
    df_gene = df_gene.transpose()

    if df_gene.shape[1] == 0:
        raise Exception("All cells have been filtered out!")

    # Apply same filtering to transcript matrix
    # Todo filter mitochondrial transcripts
    df_tr["n_transcripts"] = df_tr.astype(bool).sum(axis=1)
    n_dropped = \
        df_tr[df_tr["n_transcripts"] < args.min_genes].shape[0]
    logger.info(f"Dropping {n_dropped} cells with < {args.min_genes} genes")
    df_tr = df_tr[df_tr["n_transcripts"] >= args.min_genes]
    df_tr = df_tr.drop(["n_transcripts"], axis=1)
    df_tr = df_tr.transpose()

    return df_gene, df_tr


def filter_genes(df_gene,  args):
    """Filter genes."""
    logger = get_named_logger('ProcMatrix')
    df_gene["n_cells"] = df_gene.astype(bool).sum(axis=1)

    # Remove genes that are present in fewer than <args.min_cells> unique cells
    n_dropped = df_gene[df_gene["n_cells"] < args.min_cells].shape[0]
    logger.info(
        f"Dropping {n_dropped} genes observed in < {args.min_cells} cells")
    df_gene = df_gene[df_gene["n_cells"] >= args.min_cells]
    df_gene = df_gene.drop(["n_cells"], axis=1)

    if df_gene.shape[0] == 0:
        raise Exception("All genes have been filtered out!")

    return df_gene


def normalize(df_gene, args):
    """Normalize."""
    logger = get_named_logger('ProcMatrix')
    df_gene = df_gene.transpose()

    # cell_count / cell_total = X / <args.norm_count>
    logger.info(f"Normalizing counts to {args.norm_count} reads per cell")
    genes = df_gene.columns != "total"
    df_gene.loc[:, genes] = df_gene.loc[:, genes] * args.norm_count
    df_gene.loc[:, genes] = df_gene.loc[:, genes].div(df_gene["total"], axis=0)
    df_gene = df_gene.drop(["total"], axis=1)
    df_gene = df_gene.transpose()

    return df_gene


def main(args):
    """Run entry point."""
    logger = get_named_logger('ProcMatrix')
    df_gene = pd.read_csv(
        args.gene_counts, sep="\t").set_index("gene")
    df_transcript = pd.read_csv(
        args.transcript_counts, sep="\t").set_index("transcript")

    df_gene, df_transcript = filter_cells(df_gene, df_transcript, args)

    df_gene = filter_genes(df_gene, args)
    df_transcript = filter_genes(df_transcript, args)

    df_gene = normalize(df_gene, args)
    df_transcript = normalize(df_transcript, args)

    df_gene = np.log10(df_gene + 1)
    df_transcript = np.log10(df_transcript + 1)

    logger.info(
        f"Processed gene matrix: {df_gene.shape[0]} "
        f"genes x {df_gene.shape[1]} cells")
    df_gene.to_csv(f"{args.sample_id}.gene_expression.processed.tsv", sep="\t")

    logger.info(
        f"Processed transcript matrix: {df_transcript.shape[0]} "
        f"transcripts x {df_transcript.shape[1]} cells")
    df_transcript.to_csv(
        f"{args.sample_id}.transcript_expression.processed.tsv", sep="\t")


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
