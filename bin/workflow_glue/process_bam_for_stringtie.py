#!/usr/bin/env python
"""Preprocess bam ready for transcript assemble with stringtie.

Toggle the reverese flag to effectively invert the alignment to put in correct
orientation for stringtie.
"""
from pysam import AlignmentFile as AlnFile

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("CalSat")

    parser.add_argument(
        "bam",
        help="bam",
        default="output.png",
    )

    parser.add_argument(
        "contig",
        help="contig to process"
    )
    return parser


def main(args):
    """Entry point.

    Up to this point in the workflow the alignments are in reverse
    orientation in relation to the mRNA forsome of the 10x kits.
    This simply reverses the orientation of the alignment.
    """
    bam_file = args.bam

    with AlnFile(bam_file, "rb") as bam, AlnFile("-", "w", template=bam) as bam_out:
        for align in bam.fetch(contig=args.contig):
            align.flag ^= 16
            bam_out.write(align)


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
