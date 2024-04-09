"""Add tags from a text file to a BAM."""

from pathlib import Path

import pandas as pd
import pysam
from .util import get_named_logger, wf_parser  # noqa: ABS101

logger = get_named_logger("TagBAMs")


def argparser():
    """Create argument parser."""
    parser = wf_parser("tag_bams")

    parser.add_argument(
        "in_bam", type=Path,
        help="BAM file for tagging")

    parser.add_argument(
        "out_bam", type=Path,
        help="Path for tagged output BAM")

    parser.add_argument(
        "tags", type=Path,
        help="Read tags TSV")

    parser.add_argument(
        "chrom",
        help="Chromosome name")

    parser.add_argument(
        "--threads", default=2, type=int,
        help="Number of threads used for BAM reading/writing.")
    return parser


def add_tags(tags_file, in_bam, out_bam, chrom, threads):
    """Add all the required tags to the BAM file."""
    logger.info("Reading tag data.")

    # read everything as str since we need to serialse to string anyway
    tags = pd.read_csv(
        tags_file, sep='\t', dtype=str, index_col="read_id").rename(
            columns={"gene": "GN", "transcript": "TR"}, copy=False)
    logger.info("Indexing tag data by read_id")
    tags = tags.to_dict(orient="index")

    logger.info("Tagging reads.")
    names = "CR CB CY UR UB UY GN TR".split()
    with pysam.AlignmentFile(
            in_bam, "rb", threads=threads) as bam_in:
        with pysam.AlignmentFile(
                out_bam, "wb", template=bam_in, threads=threads) as bam_out:
            for align in bam_in.fetch(contig=chrom):
                read_id = align.query_name
                try:
                    row = tags[read_id]
                except KeyError:
                    continue  # don't write reads without tags
                else:
                    for key in names:
                        align.set_tag(key, row[key], value_type="Z")
                    bam_out.write(align)
    logger.info("Done.")


def main(args):
    """Entry point."""
    add_tags(args.tags, args.in_bam, args.out_bam, args.chrom, args.threads)
