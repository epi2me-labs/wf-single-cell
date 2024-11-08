"""Prepare data for the report."""
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.io import mmread

from .adapter_scan_vsearch import AdapterSummary  # noqa: ABS101
from .create_matrix import ExpressionSummary  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


# TODO: Code in this script move into the main report.py script


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("prepare_report_data")

    parser.add_argument(
        "sample_id",
        help="ID of the sample being processed")
    parser.add_argument(
        "adapter_stats", type=Path,
        help="Workflow summary statistics")
    parser.add_argument(
        "bam_stats", type=Path,
        help="Alignment summary statistics")
    parser.add_argument(
        "expression_stats", type=Path,
        help="Expression summary statistics")
    parser.add_argument(
        "white_list",
        help="Workflow summary statistics")
    parser.add_argument(
        "survival_out", type=Path,
        help="Output TSV with survival data for each stage.")
    parser.add_argument(
        "bam_stats_out", type=Path,
        help="Output TSV with combined alignment summary stats.")
    parser.add_argument(
        "raw_gene_expression", type=Path,
        help="Sparse data in MEX format.")
    parser.add_argument(
        "genes_of_interest", type=Path,
        help="TSV file of file names.")
    parser.add_argument(
        "n_input_seqs", type=int,
        help="Number of seqs input to the workflow after read quality filtering.")
    return parser


def combine_bam_stats(input_dir, sample_id):
    """Aggregate alignment statistics."""
    dfs = []
    colnames = {
        "PrimAln": "primary",
        "SecAln": "secondary",
        "SupAln": "supplementary",
        "Unmapped": "unmapped",
        "TotalReads": "reads_aligned"

    }
    for stats in input_dir.glob('*.tsv'):
        dfs.append(pd.read_csv(
            stats, sep='\t',
            usecols=colnames.keys(),
            dtype=int
        ))
    df = pd.concat(dfs)
    df = pd.DataFrame(df.sum(axis=0)).T
    df = df.rename(columns=colnames)
    df.insert(0, 'sample', sample_id)
    df.insert(1, 'reads_aligned', df.pop('reads_aligned'))

    return df


def combine_expression_stats(input_dir):
    """Summarise expressions summary files."""
    fnames = list(input_dir.glob("*.json"))
    if len(fnames) == 0:
        raise IOError("No summary JSON files found.")

    summary = ExpressionSummary.from_json(fnames[0])
    if len(fnames) > 1:
        for other in fnames[1:]:
            summary += ExpressionSummary.from_json(other)
    return summary


def combine_adapter_stats(input_dir):
    """Combine adapter configuration summary files."""
    fnames = list(input_dir.glob("*.json"))
    if len(fnames) == 0:
        raise IOError("No summary JSON files found.")

    summary = AdapterSummary.from_json(fnames[0])
    if len(fnames) > 1:
        for other in fnames[1:]:
            summary += AdapterSummary.from_json(other)
    return summary


def get_total_cells(white_list):
    """Create dataframe with total cells."""
    # ok this is a little cheesy, but consistent for ease
    total_cells = len(pd.read_csv(white_list, sep='\t', header=None))
    return {"cells": total_cells}


def get_genes_of_interest_expression(mex_dir, genes):
    """Get a subset of the expression data.

    Given a list of genes, extract corresponding expression data from the MEX format
    matrix.
    """
    genes_to_plot = pd.read_csv(genes, header=None)[0]
    matrix = mmread(mex_dir / 'matrix.mtx.gz')
    barcodes = pd.read_csv(mex_dir / 'barcodes.tsv.gz', header=None)
    # Remove '-1' suffix from barcodes
    barcodes = barcodes[0].str.split('-', expand=True)[0]
    features = pd.read_csv(mex_dir / 'features.tsv.gz', sep='\t', header=None)[1]
    rows = []
    for g in genes_to_plot:
        try:
            feature_idx = features[features == g].index[0]
            rows.append(np.array(matrix.getrow(feature_idx).todense()).flatten())
        except IndexError:
            rows.append([np.nan] * len(barcodes))

    return (
        pd.DataFrame.from_records(rows, index=genes_to_plot, columns=barcodes)
    )


def main(args):
    """Entry point for script."""
    logger = get_named_logger('PrepReport')
    logger.info('Preparing report data.')
    stats = dict()
    stats.update(combine_expression_stats(args.expression_stats))
    stats.update(combine_adapter_stats(args.adapter_stats))
    stats.update(get_total_cells(args.white_list))
    # n seqs after any read quality filtering
    n_input_reads = args.n_input_seqs
    stats.update({'reads': n_input_reads})

    survival = (
        pd.DataFrame.from_dict(stats, orient="index", columns=['count'])
        .reset_index(names="statistic"))

    # this is a little nonsensical for some stats
    survival['Proportion of reads [%]'] = 100 * survival['count'] / n_input_reads
    survival['sample_id'] = args.sample_id
    survival.to_csv(args.survival_out, sep='\t', index=False)

    aln_stats = combine_bam_stats(args.bam_stats, args.sample_id)
    aln_stats.to_csv(args.bam_stats_out, sep='\t', index=False)

    goi_df = get_genes_of_interest_expression(
        args.raw_gene_expression, args.genes_of_interest)
    goi_df.to_csv(
        Path(f"{args.sample_id}_expression") / 'raw_goi_expression.tsv', sep='\t')
