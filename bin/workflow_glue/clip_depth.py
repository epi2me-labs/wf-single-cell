"""Clip read depth."""

import random

import pandas as pd
from pysam import AlignmentFile
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("clip_depth")

    parser.add_argument(
        "--bed",
        help="Regions to clip")

    parser.add_argument(
        "--bam_in",
        help="input bam file")

    parser.add_argument(
        "--target_depth",
        help="Desired read depth",
        type=int)

    parser.add_argument(
        "--bam_out",
        help="output bam file")

    return parser


def main(args):
    """Run entry point."""
    df_hi_cov = pd.read_csv(
        args.bed, sep='\t',
        names=['read_id', 'depth'],
        index_col='read_id',
        dtype={
            'read_id': str,
            'depth': float
        }
    )

    random.seed(1889)
    with AlignmentFile(args.bam_in, "rb", check_sq=False) as bam:
        with AlignmentFile(args.bam_out, "wb", template=bam) as out_bam:
            for aln in bam.fetch(until_eof=True):
                if aln.query_name in df_hi_cov.index:
                    window_depth = df_hi_cov.at[aln.query_name, 'depth']
                    # Randomly discard records to achieve target depth
                    # Random returns a random floating number between 0 and 1
                    if random.random() > args.target_depth / window_depth:
                        continue  # Discard record
                out_bam.write(aln)
