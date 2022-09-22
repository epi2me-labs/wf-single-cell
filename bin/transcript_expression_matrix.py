#!/usr/bin/env python
"""Build transcript expression matrix.

Take .tmap file from stringtie (concatenated from all barcode files)
and rejig to transcript x cell matrix of TPM values.

Filter assembled transcripts based on assigned gffcompare classes.
     i: Fully contained in a reference intron
     p: Possible polymerase run-on
     u: novel transcripts
     s: intron match on opposite strand

Filter genes and barcodes based on what is present in the gene matrix, which
has already been filtered for lowly-expressed genes and rare barcodes.

Novel transcripts are removed as they are currently generated per barcode
and cannot be compared with each other.

TODO: group togeter novel transcripts from each barcode and give the same id
so that tehy can be comapred across the barcodes.
"""

import argparse
import logging

import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    """Create argument parser."""
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "tmap",
        help="Concatenated results of gffcomapre .tmap file.",
        type=str,
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
        "--norm_count",
        help="Normalize to this number of counts per cell [10000]",
        type=int,
        default=10000,
    )

    parser.add_argument(
        "--output",
        help="Output TSV file containing processed gene expression matrix, \
               where genes are rows and cells are columns \
               [expression.processed.tsv]",
        type=str,
        default="expression.processed.tsv",
    )

    parser.add_argument(
        "--gene_matrix",
        help="The already processed gene matrix. This will be used to filter \
             the cells and transcripts",
        type=str
    )

    parser.add_argument(
        "--mito_output",
        help="Output TSV file containing percentage of UMI counts coming from \
               mitochondrial genes [expression.mito.tsv]",
        type=str,
        default="expression.mito.tsv",
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: 2=info, 3=warning [2]",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    """Initiate logger."""
    logging.basicConfig(
        format="%(asctime)s -- %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def main(args):
    """Entry point."""
    init_logger(args)

    header = [
        'ref_gene_id',
        'ref_id',
        'class_code',
        'qry_gene_id',
        'qry_id',
        'num_exons',
        'FPKM',
        'TPM',
        'cov',
        'len',
        'major_iso_id',
        'ref_match_len',
        'barcode']

    df_tmap = pd.read_csv(args.tmap, sep='\t', index_col=None, header=None)
    df_tmap.columns = header

    gene_df = pd.read_csv(args.gene_matrix, sep='\t')
    df_tmap = filter_transcripts(gene_df, df_tmap)
    df_tmap = filter_cells(gene_df, df_tmap)
    df_tmap = filter_gffcompare_classes(df_tmap)

    rows = []

    grouped = df_tmap.groupby('ref_id')

    # Create transcript x cell matrix for each barcode, then concatenate.
    for _, dfg in grouped:
        # Barcode/transcript duplicates can exist if multiple query
        # transctipts made by stringtie map to the same reference transcript
        # If there are any duplicates, keep one and replace the TPM count with
        # The summed from all duplicates.
        dfg = dfg.reset_index(drop=True)
        dups = dfg.duplicated(['barcode', 'ref_id'], keep=False)

        if any(dups):
            dup_df = dfg.loc[dups]
            summed = dup_df.groupby(['barcode', 'ref_id']).sum().reset_index()
            dupes_to_remove = dfg.duplicated(['barcode', 'ref_id'])
            dfg = dfg[~dupes_to_remove]
            # Put the summed values back
            for _, row in summed.iterrows():
                dfg.loc[
                    (dfg.barcode == row.barcode) &
                    (dfg.ref_id == row.ref_id), 'TPM'] = row.TPM

        dfp = dfg.pivot(index='ref_id', columns='barcode', values='TPM')
        rows.append(dfp)

    result_df = pd.concat(rows)

    result_df.to_csv(args.output, sep='\t')


def filter_gffcompare_classes(df):
    """Filter some gffcompare transcript classes."""
    df = df[~df['class_code'].isin(['i', 'p', 'u', 's'])]
    return df


def filter_cells(df_genes, df_tmap):
    """Remove cells not in gene matrix."""
    cells = df_genes.columns
    df_tmap = df_tmap[df_tmap['barcode'].isin(cells)]
    return df_tmap


def filter_transcripts(df_genes, df_tmap):
    """Remove transcripts whose parent gene is not in the gene matrix."""
    genes = df_genes['gene']
    df_tmap = df_tmap[df_tmap['ref_gene_id'].isin(genes)]
    return df_tmap


if __name__ == "__main__":
    args = parse_args()

    main(args)
