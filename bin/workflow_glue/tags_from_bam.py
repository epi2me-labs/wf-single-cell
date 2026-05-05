
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
                        # SpaceRanger 4.x Visium HD writes the raw 28-mer
                        # sequence barcode in CB and the spatial-form barcode
                        # (s_002um_NNNNN_NNNNN-1) in the `sb` tag. The
                        # wf-single-cell HD codepath (is_visium_hd check,
                        # 2um->8um binning, spatial plots in the report) keys
                        # off the spatial form, so substitute `sb` for `CB`
                        # when present. Without this the HD branch is
                        # silently bypassed and `makeReport` is skipped
                        # because gene_expression_8um is empty.
                        if record.has_tag('sb'):
                            val = record.get_tag('sb')
                        barcode_counter[val] += 1
                    # Sanitize string tag values so the TSV is parseable by
                    # pandas (default quotechar='"'). Quality strings (CY/UY)
                    # routinely contain ASCII 0x22 ('"' = Phred Q1) which
                    # otherwise causes "EOF inside string" errors in
                    # downstream pd.read_csv. Also strip tab/CR/LF which
                    # would corrupt TSV row/column structure.
                    if isinstance(val, str):
                        val = (val.replace('"', "'")
                                  .replace('\t', ' ')
                                  .replace('\r', ' ')
                                  .replace('\n', ' '))
                else:
                    val = '-'
                tag_values.append(val)

            tags_str = f"{record.query_name}\t" + "\t".join(tag_values)
            tags_out.write(f"{tags_str}\n")

    df_bc_counts = (
        pd.DataFrame(
            barcode_counter.items(), columns=['barcode', 'count'])
        .sort_values(by='count', ascending=False)
    )

    if len(df_bc_counts) > 0:
        # Bin barcode counts to 8um coordinates

        # Extract coordinate part from barcode. SpaceRanger BAMs can carry CB
        # tags that don't follow the Visium HD coordinate pattern
        # `_NNNNN_NNNNN-` (e.g. on decoy/unplaced contigs). Drop those rows
        # rather than crashing on NaN -> int conversion. They cannot be placed
        # spatially by definition, so excluding them from the binning summary
        # has no effect on per-read or per-gene downstream output.
        coords = df_bc_counts['barcode'].str.extract(r'_(\d{5})_(\d{5})-')
        valid = coords.notna().all(axis=1)
        n_dropped = int((~valid).sum())
        if n_dropped:
            logger.warning(
                f"Dropped {n_dropped} barcode(s) from spatial binning that "
                f"don't match the Visium HD coordinate pattern "
                f"_NNNNN_NNNNN- (only affects barcode_counts.tsv; per-read "
                f"tags TSV is unaffected).")
        df_bc_counts = df_bc_counts.loc[valid].copy()
        df_bc_counts[["x_2um", "y_2um"]] = coords.loc[valid].astype(int)

        # Convert to 8um coordinates by integer division (bin size = 4)
        df_bc_counts["x_8um"] = df_bc_counts["x_2um"] // 4
        df_bc_counts["y_8um"] = df_bc_counts["y_2um"] // 4

        df_bc_counts.drop(columns=["x_2um", "y_2um"], inplace=True)

        # Group by 8um bins and sum counts
        df_bc_counts = \
            df_bc_counts.groupby(["x_8um", "y_8um"])["count"].sum().reset_index()

        # Generate a new 8um bin barcode name eg s_008_um_01234_01234
        df_bc_counts["barcode"] = df_bc_counts.apply(
            lambda row: (
                f's_008_um_{str(row["x_8um"]).zfill(5)}_{str(row["y_8um"]).zfill(5)}'),
            axis=1
        )

        df_bc_counts.drop(columns=["x_8um", "y_8um"], inplace=True)
        df_bc_counts = df_bc_counts[["barcode", "count"]]

    # Output to TSV
    df_bc_counts.to_csv(
        args.barcode_counts_out, header=True, sep="\t", index=False)
