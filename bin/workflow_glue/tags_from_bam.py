
"""Extract tags from BAM file and write to TSV."""
from collections import Counter

import pandas as pd
import pysam
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("prepare_report_data")

    parser.add_argument(
        'bam_in', help="Input FASTQ file. Can be path or stdin (-)")
    parser.add_argument(
        'tags_out', help="Output TSV file with tags per read")
    parser.add_argument(
        'barcode_counts_out', help="Output TSV file with tags per read")
    parser.add_argument(
        "--tags", nargs='*',
        help="List of tags to extract from the BAM file")
    parser.add_argument(
        "--chrom", default=None,
        help="List of tags to extract from the BAM file")

    return parser


def main(args):
    """Read in BAM and extract tags. Write TSV of tags per read.."""
    logger = get_named_logger("TagsFromBAM")
    logger.info("Extracting tags from BAM file")

    # Write output file header
    with open(args.tags_out, "w") as tags_out:
        tags_out.write(
            "read_id\t" + "\t".join(args.tags) + "\n")

    barcode_counter = Counter()

    in_ = pysam.AlignmentFile(args.bam_in, "rb")
    out = open(args.tags_out, "a")

    with in_ as bam_in, out as tags_out:
        if args.chrom:
            iterator = bam_in.fetch(args.chrom, until_eof=True)
        else:
            iterator = bam_in.fetch(until_eof=True)
        for record in iterator:
            tag_values = []

            if not record.has_tag('CB'):
                continue

            for tag in args.tags:
                if record.has_tag(tag):
                    val = record.get_tag(tag)
                    if tag == 'CB':
                        barcode_counter[val] += 1
                    # Quote double quaotes in quality strings
                    # For now just report ?
                    if tag in ['CY', 'UY']:
                        # val = val.replace('"', '""')
                        # val = f'"{val}"'
                        val = "?" * 8
                else:
                    val = '-'
                tag_values.append(val)

            tags_str = f"{record.query_name}\t" + "\t".join(tag_values)
            tags_out.write(f"{tags_str}\n")

    (
        pd.DataFrame(
            barcode_counter.items(), columns=['barcode', 'count'])
        .sort_values(by='count', ascending=False)
        .to_csv(args.barcode_counts_out, sep='\t', index=False, header=True)
    )
