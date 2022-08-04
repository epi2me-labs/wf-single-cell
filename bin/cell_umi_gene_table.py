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
        "bam",
        help="Aligned BAM file with gene, barcode, and UMI tags",
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


def load_bam_entries(n_aligns, args):
    """
    Load read_id, gene, corrected barcode, and corrected UMI values from BAM
    """
    bam = pysam.AlignmentFile(args.bam, "rb")

    records = []
    for i, align in enumerate(tqdm(bam.fetch(), total=n_aligns)):
        read_id = align.query_name

        # Annotated gene name = GN:Z
        gene = align.get_tag("GN")
        # Corrected cell barcode = CB:Z
        bc = align.get_tag("CB")
        # Corrected UMI = UB:Z
        umi = align.get_tag("UB")

        records.append((read_id, gene, bc, umi))

    bam.close()

    df = pd.DataFrame.from_records(
        records, columns=["read_id", "gene", "barcode", "umi"]
    )

    return df


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
    n_aligns = int(sum([contig.mapped for contig in stats]))
    chroms = dict(
        [(contig.contig, contig.mapped) for contig in stats if contig.mapped > 0]
    )
    bam.close()
    return n_aligns, chroms


def main(args):
    init_logger(args)

    logger.info(f"Counting alignments in {args.bam}")
    n_aligns, chroms = get_bam_info(args.bam)

    logger.info(f"Reading read_ids, genes, barcodes, and UMIs from {args.bam}")
    df = load_bam_entries(n_aligns, args)

    logger.info(f"Writing data to {args.output}")
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    args = parse_args()

    main(args)
