#!/usr/bin/env python
"""Add gene tags."""
import os
from pathlib import Path

import numpy as np
import pysam
from tqdm import tqdm

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("add_gene_tags")

    # Positional mandatory arguments
    parser.add_argument("bam", help="Sorted BAM file", type=Path)

    parser.add_argument(
        "gene_assigns",
        help="TSV read/gene assignments file. \
        IMPORTANT: reads in the input BAM and gene_assigns file must have the \
        same order.",
        type=Path,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output BAM file containing aligned reads with gene name \
            tags (GN) [gene.sorted.bam]",
        type=Path,
        default="gene.sorted.bam",
    )

    return parser


def get_bam_info(bam):
    """Get number of alignments and names of all contigs in the reference.

    :param bam: Path to sorted BAM file
    :type bame: str
    :return: Sum of all alignments in the BAM index file and list of all chroms
    :rtype: int,list
    """
    with pysam.AlignmentFile(bam, "rb") as bam:
        stats = bam.get_index_statistics()
        n_aligns = int(np.sum([contig.mapped for contig in stats]))
        chroms = [contig.contig for contig in stats]

    return n_aligns, chroms


def process_bam_entries(args):
    """Process bam entries.

    Iterate simultaneously throught the BAM and the featureCounts TSV file of
    read/gene assignments. These files must be sorted identically, otherwise it
    will throw an exception. The gene names found in the featureCounts TSV are
    added to the GN tag for each alignment in the output BAM.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    logger = get_named_logger('AddGeneTag')
    n_reads, chroms = get_bam_info(args.bam)

    with pysam.AlignmentFile(args.bam, "rb") as bam:
        with pysam.AlignmentFile(args.output, "wb", template=bam) as bam_out:

            # If input BAM file is empty or there are no gene assignments,
            # write an empty output BAM file
            if (n_reads == 0) or (os.path.getsize(args.gene_assigns) == 0):
                pass
            else:
                logger.info(f"Adding gene tags (GN) to {args.bam}")
                with open(args.gene_assigns, "r") as f:
                    for align in tqdm(bam.fetch(), total=n_reads):
                        line = f.readline().strip()
                        gene_assigns_id = line.split("\t")[0]
                        gene_assigns_gene = line.split("\t")[3]

                        assert (
                            gene_assigns_id == align.query_name
                        ), "BAM and featureCounts reads not ordered"

                        # Annotated gene name = GN:Z
                        align.set_tag("GN", gene_assigns_gene, value_type="Z")

                        bam_out.write(align)


def main(args):
    """Run the entry point."""
    process_bam_entries(args)


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
