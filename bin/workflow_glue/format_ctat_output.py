"""Convert ctat-LR-fusion outputs from long to short format."""
from collections import defaultdict
import csv

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Parse the arguments."""
    parser = wf_parser(
        "Map fusions to cell barcodes.")
    parser.add_argument(
        "fusion_file",
        help="Path to the fusion output file (TSV format).")
    parser.add_argument(
        "read_info_file",
        help="Path to the per-read info file (TSV format).")
    parser.add_argument(
        "per_read_output",
        help="Path to save the output TSV file.")
    parser.add_argument(
        "per_fusion_output",
        help="Path to save the output TSV file.")
    parser.add_argument(
        "cell_summary_out",
        help="Path to save the output TSV file.")
    parser.add_argument(
        "sample_id",
        help="sample identifier to add to tables.")
    parser.add_argument(
        "--unmatched_reads_out",
        help="Path to save the output TSV file of fusions with no matching CB.",
        default='unmatched_reads.txt')

    return parser


def load_fusion_data(fusion_file):
    """Load fusion data. Extract relevant fields along with read associations."""
    logger = get_named_logger('FmtCtat')
    try:
        per_fusion_df = pd.read_csv(fusion_file, sep="\t")
    except pd.errors.EmptyDataError:
        logger.warning(
            f"""The fusion file {fusion_file} is empty.
            No candidate fusions found by ctat-LR-fusion.""")
        return None
    if len(per_fusion_df) == 0:
        logger.warning(
            f"""The fusion file {fusion_file} contained no entries.
            No candidate fusions passed ctat-LR-fusion filters.""")
        return None
    per_fusion_df.rename(columns={"#FusionName": "FusionName"}, inplace=True)

    # Convert read IDs to lists
    per_fusion_df["LR_accessions"] = per_fusion_df["LR_accessions"].str.split(",")

    # Expand multiple read IDs per fusion
    per_read_df = per_fusion_df.explode("LR_accessions")
    # Select relevant columns
    per_read_df = per_read_df[[
        "FusionName", "LeftGene", "LeftBreakpoint", "RightGene",
        "RightBreakpoint", "SpliceType", "LR_accessions"
    ]].rename(columns={"LR_accessions": "read_id"})

    # Handle duplicate read IDs using defaultdict (list)
    fusion_dict = defaultdict(list)
    for row in per_read_df.itertuples():
        fusion_dict[row.read_id].append({
            "FusionName": row.FusionName,
            "LeftGene": row.LeftGene,
            "LeftBreakpoint": row.LeftBreakpoint,
            "RightGene": row.RightGene,
            "RightBreakpoint": row.RightBreakpoint,
            "SpliceType": row.SpliceType
        })

    logger.info(f"Total fusions processed: {per_fusion_df['FusionName'].nunique()}")
    logger.info(f"Total unique reads linked to fusions: {len(fusion_dict)}")

    return fusion_dict


def process_read_info(read_info_file, fusion_dict):
    """Combine single-cell tags with fusion info."""
    matched_results = []
    unmatched_reads = set(fusion_dict.keys())  # Track reads missing barcode/UMI

    with open(read_info_file, newline='') as csvfile:
        reader = csv.DictReader(csvfile, delimiter="\t")

        for line in reader:
            read_id = line['read_id']
            if read_id in fusion_dict:  # Only process relevant reads
                for fusion in fusion_dict[read_id]:
                    matched_results.append({
                        "FusionName": fusion["FusionName"],
                        "LeftGene": fusion["LeftGene"],
                        "LeftBreakpoint": fusion["LeftBreakpoint"],
                        "RightGene": fusion["RightGene"],
                        "RightBreakpoint": fusion["RightBreakpoint"],
                        "SpliceType": fusion["SpliceType"],
                        "CB": line["corrected_barcode"],
                        "UB": line["corrected_umi"],
                        "read_id": read_id
                    })
                unmatched_reads.discard(read_id)  # Remove matched reads

    return pd.DataFrame(matched_results), unmatched_reads


def main(args):
    """Run the script."""
    logger = get_named_logger('FmtCtat')

    logger.info("Loading fusion data...")

    fusion_dict = load_fusion_data(args.fusion_file)

    if fusion_dict is None:
        logger.warning("No fusion data found. writing empty files.")

        with open(args.per_read_output, 'w') as fh1:
            fh1.write(
                "FusionName\tLeftGene\tLeftBreakpoint\tRightGene\tRightBreakpoint"
                "\tSpliceType\tCB\tUB\tread_id\n"
            )

        with open(args.per_fusion_output, 'w') as fh2:
            fh2.write(
                "Fusion\tLeftGene\tLeftBreakpoint\tRightGene\tRightBreakpoint"
                "\tSpliceType\tcells\tUMIs\tsample_ID\n"
            )
        (
            pd.DataFrame.from_records(
                [[args.sample_id, 0, 0, 0, 0, 0]],
                columns=[
                    'sample_ID', 'cells_with_fusions', 'unique_fusions', 'reads',
                    'mean_fusion_reads_per_cell', 'mean_unique_fusions_per_cell'],
            )
            .to_csv(args.cell_summary_out, sep="\t", index=False)
        )
    else:
        logger.info("Processing read information...")
        merged_df, unmatched_reads = process_read_info(
            args.read_info_file, fusion_dict)

        # Make per-fusion summary from barcode-assigned fusion + reads.
        (
            merged_df.groupby(
                ['FusionName', 'LeftGene', 'LeftBreakpoint',
                 'RightGene', 'RightBreakpoint', 'SpliceType'])
            .agg(
                cells=('CB', 'nunique'),
                UMIs=('UB', 'nunique'))
            .reset_index()
            .assign(
                sample_ID=args.sample_id)
            .rename(columns={'FusionName': 'Fusion'})
            .to_csv(args.per_fusion_output, sep="\t", index=False)
        )

        # Sort reads by most commonly occuring fusion pair
        merged_df = (merged_df.sort_values(
                by="FusionName",
                key=lambda col: col.map(col.value_counts()), ascending=False))

        logger.info("Saving merged single cell/fusion output")
        merged_df.to_csv(args.per_read_output, sep="\t", index=False)

        logger.info("Writing the per fusion summary...")
        # Regerate the per-fusion summary adding some cell-specific info.
        (
            pd.DataFrame.from_dict({
                'sample_ID': args.sample_id,
                'cells_with_fusions': merged_df['CB'].nunique(),
                'unique_fusions': merged_df['FusionName'].nunique(),
                'reads': len(merged_df),
                'mean_fusion_reads_per_cell': (
                    merged_df.groupby('CB')['FusionName'].count().mean()),
                'mean_unique_fusions_per_cell': (
                    merged_df.groupby('CB')['FusionName'].nunique().mean()),
                }, orient='index')
            .T
            .to_csv(args.cell_summary_out, sep="\t", index=False)
        )

        # Summary of unmatched reads
        logger.info(
            (
                "\n**Summary:**"
                f"Matched reads with barcode/UMI: {len(merged_df)}\n"
                f"Reads missing barcode/UMI info: {len(unmatched_reads)}")
            )

        if unmatched_reads:
            with open(args.unmatched_reads_out, "w") as uf:
                for read in unmatched_reads:
                    uf.write(f"{read}\n")
        logger.info("Process complete!")
