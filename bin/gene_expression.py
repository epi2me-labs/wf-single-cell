import argparse
import collections
import gzip
import logging
import os
import pathlib
import re
import subprocess
import sys

import numpy as np
import pandas as pd
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument("bam", help="Sorted BAM file", type=str)

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output TSV file containing gene expression matrix, where genes \
        are rows and barcodes are columns [gene_expression.tsv]",
        type=str,
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
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def get_bam_info(bam):
    """
    Use `samtools idxstat` to get number of alignments and names of all contigs
    in the reference.

    :param bam: Path to sorted BAM file
    :type bame: str
    :return: Sum of all alignments in the BAM index file and list of all chroms
    :rtype: int,list
    """
    bam = pysam.AlignmentFile(bam, "rb")
    stats = bam.get_index_statistics()
    n_aligns = int(np.sum([contig.mapped for contig in stats]))
    chroms = [contig.contig for contig in stats]
    bam.close()
    return n_aligns, chroms


def read_bam_entries(bam, n_reads):
    """ """
    genes = set()
    barcodes = set()

    umi_sets = collections.defaultdict(set)

    # Build regular expression for "gene" annotations where no gene was found
    REGEX = r"[a-zA-Z0-9]+_\d+_\d+"  # e.g. chr7_44468000_44469000
    for align in tqdm(bam.fetch(), total=n_reads):
        gene = align.get_tag("GN")
        barcode = align.get_tag("CB")
        umi = align.get_tag("UB")

        m = re.search(REGEX, gene)
        if m is None:
            genes.add(gene)
            umi_sets[(gene, barcode)].add(umi)

        barcodes.add(barcode)

    return genes, barcodes, umi_sets


def populate_matrix(genes, barcodes, umi_sets):
    """ """
    rows = list(genes)
    rows.sort()

    cols = list(barcodes)
    cols.sort()

    df = pd.DataFrame(0, index=rows, columns=cols)

    for (gene, barcode), umis in tqdm(umi_sets.items()):
        df.loc[gene, barcode] = len(umis)

    return df


def process_bam_entries(args):
    """
    Iterate through the BAM file and count unique UMIs (UB tag) associated with
    each gene (GN tag) and cell barcode (CB tag).

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    n_reads, chroms = get_bam_info(args.bam)

    bam = pysam.AlignmentFile(args.bam, "rb")

    # If input BAM file is empty or there are no gene assignments, raise exception
    assert n_reads > 0, "WARNING: no alignments detected!"

    logger.info(f"Building gene expression matrix from {args.bam}")
    genes, barcodes, umi_sets = read_bam_entries(bam, n_reads)

    logger.info(f"Populating gene expression matrix from {args.bam}")
    df = populate_matrix(genes, barcodes, umi_sets)

    df.to_csv(args.output, sep="\t", index_label="gene")


def main(args):
    init_logger(args)

    process_bam_entries(args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
