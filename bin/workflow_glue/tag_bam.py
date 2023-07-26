#!/usr/bin/env python
"""Cluster UMIs."""

from pathlib import Path

import pandas as pd
import pysam
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("tag_bams")

    parser.add_argument(
        "--in_bam",
        help="BAM file for tagging",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--out_bam",
        help="Path for tagged output BAM",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--tags",
        help="Read tags TSV",
        type=Path,
        required=True
    )

    parser.add_argument(
        "--chrom",
        help="Chromosome name",
        required=True
    )
    parser.add_argument(
        "--flip",
        help="Reverse the seq by toggling the SAM flag",
        action='store_true'
    )
    return parser


def add_tags(tags_file, in_bam, out_bam, chrom, rev):
    """Add all the required tags to the BAM file."""
    df = pd.read_csv(tags_file, sep='\t', index_col=0)
    with pysam.AlignmentFile(in_bam, "rb") as bam_in:
        with pysam.AlignmentFile(out_bam, "wb", template=bam_in) as bam_out:
            for align in bam_in.fetch(contig=chrom):
                read_id = align.query_name
                try:
                    row = df.loc[read_id]
                except KeyError:
                    continue  # No barcode/umi for this read
                # uncorrectred cell barcode
                align.set_tag('CR', row['CR'], value_type="Z")
                # Correctred cell barcode
                align.set_tag('CB', row['CB'], value_type="Z")
                # barcode qscores
                align.set_tag('CY', row['CY'], value_type="Z")
                # Uncorrected UMI
                align.set_tag('UR', row['UR'], value_type="Z")
                # UMI quality score
                align.set_tag('UY', row['UY'], value_type="Z")
                # Corrected UMI = UB:Z
                align.set_tag("UB", row['UB'], value_type="Z")
                # Annotated gene name = GN:Z
                align.set_tag("GN", row['gene'], value_type="Z")
                # Annotated transcript name = TR:Z
                align.set_tag("TR", row['transcript'], value_type="Z")

                # Up to this point in the workflow the alignments are in reverse
                # orientation in relation to the mRNA for 31 and multiome kits.
                # Optionally reverse the orientation for the output BAMs.
                if rev:
                    align.flag ^= 16  # reverse read alignment flag

                bam_out.write(align)


def main(args):
    """Entry point."""
    add_tags(args.tags, args.in_bam, args.out_bam, args.chrom, args.flip)
