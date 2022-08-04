import argparse
import gzip
import logging
import math
import multiprocessing
import os
import pathlib
import re
import shutil
import sys
import tempfile

import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing reads aligned to a reference.",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "--output_dir",
        help="Output directory where chromosome-specific BAM files will be \
        written. [./chrom_bams]",
        type=str,
        default="./chrom_bams",
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


def process_bam_records(tup):
    """
    Write a separate BAM file for alignments to each chromosome in the
    input_bam.

    :param tup: Tuple containing the input arguments
    :type tup: tup
    :return: Path to a temporary BAM file
    :rtype: str
    """
    input_bam = tup[0]
    chrom = tup[1]
    args = tup[2]

    # Open input BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Open output BAM file for writing
    chrom_bam_fn = os.path.join(args.output_dir, f"{chrom}.sorted.bam")
    bam_out = pysam.AlignmentFile(chrom_bam_fn, "wb", template=bam)

    for align in bam.fetch(contig=chrom):
        bam_out.write(align)

    bam.close()
    bam_out.close()

    # Index each newly created chrom-specific BAM
    pysam.index(chrom_bam_fn)

    return chrom_bam_fn


def launch_pool(func, func_args, procs=1):
    """
    Use multiprocessing library to create pool and map function calls to
    that pool

    :param procs: Number of processes to use for pool
    :type procs: int, optional
    :param func: Function to exececute in the pool
    :type func: function
    :param func_args: List containing arguments for each call to function <funct>
    :type func_args: list
    :return: List of results returned by each call to function <funct>
    :rtype: list
    """
    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


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
    # logger.info("Getting BAM statistics")
    n_reads, chroms = get_bam_info(args.bam)

    # Create output directory
    if os.path.exists(args.output_dir):
        shutil.rmtree(args.output_dir, ignore_errors=True)
    os.mkdir(args.output_dir)

    # Write separately BAM file for each chromosome in args.output_dir
    logger.info(f"Writing chrom-specfic BAM files to {args.output_dir}")
    func_args = []
    chroms_sorted = dict(sorted(chroms.items(), key=lambda item: item[1]))
    for chrom in chroms_sorted.keys():
        func_args.append((args.bam, chrom, args))

    launch_pool(process_bam_records, func_args, args.threads)


if __name__ == "__main__":
    args = parse_args()

    main(args)
