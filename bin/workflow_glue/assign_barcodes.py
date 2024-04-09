"""Assign barcodes.

Given a whitelist of barcodes assign raw barcodes to nearest match.
"""
import collections
from pathlib import Path

import pandas as pd
import rapidfuzz
from rapidfuzz.process import extract

from .util import get_named_logger, wf_parser  # noqa: ABS101

logger = get_named_logger("AsgnBrcdes")


def argparser():
    """Create argument parser."""
    parser = wf_parser("assign_barcodes")

    parser.add_argument(
        "whitelist", type=Path,
        help="File containing list of expected cell barcodes.")

    parser.add_argument(
        "barcode_tags", type=Path,
        help="TSV file of read_id, uncorrected_barcode, qscores.")

    parser.add_argument(
        "output_tags", type=Path,
        help="Output TSV containing columns from `barcode_tags` \
            and additional a CB (corrected barcode) column.")

    parser.add_argument(
        "output_counts", type=Path,
        help="Output TSV file containing counts for each of the assigned \
            barcodes.")

    parser.add_argument(
        "--chunksize", type=int, default=100000,
        help="Process the BAM in chunks no larger than this.")

    parser.add_argument(
        "--max_ed", type=int, default=2,
        help="Max. edit distance between putative barcode \
            and the matching whitelist barcode.")

    parser.add_argument(
        "--min_ed_diff", type=int, default=2,
        help="Min. difference in edit distance between the \
            best and second best whitelist matches.")

    return parser


def determine_barcode(bc_uncorr, whitelist, max_ed, min_ed_diff):
    """Find barcode in a whitelist corresponding to read barcode.

    :param bc_uncorr: uncorrected barcode.
    :param whitelist: list of possible barcodes.
    :param max_ed: max. edit distance between barcode and whitelist hit.
    :param min_ed_diff: min. edit distance difference between first and
        second best hits in order to accept the first as valid.
    """
    # quick return
    if bc_uncorr in whitelist:
        return bc_uncorr

    # now do levenstein on anything left
    # result is a list of tuples (bc, ed, idx) sorted by ed.
    result = extract(
        bc_uncorr,
        whitelist,
        scorer=rapidfuzz.distance.Levenshtein.distance,
        score_cutoff=max_ed + min_ed_diff + 1)

    if len(result) > 0:
        bc_match = result[0][0]
        bc_match_ed = result[0][1]
    else:
        bc_match = "X" * len(bc_uncorr)
        bc_match_ed = len(bc_uncorr)
    if len(result) > 1:
        next_match_diff = result[1][1] - bc_match_ed
    else:
        next_match_diff = len(bc_uncorr)

    # are we better than the second place?
    corrected = "-"
    if (bc_match_ed <= max_ed) and (next_match_diff >= min_ed_diff):
        corrected = bc_match
    return corrected


def process_records(
        barcode_tags, whitelist, max_ed, min_ed_diff, tags_output,
        chunksize=100000):
    """Process read barcodes stored in text file to find whitelist equivalents.

    :param barcode_tags: path to TSV with tag data
    :param whitelist: list of potential barcodes.
    :param: max_ed: max allowed edit distance between am uncorrected barcode
        and a potential corected whitelist barcode.
    :param: min_ed_diff: minimum allowed edit distance between top two
        barcode candidates.
    """
    barcode_counter = collections.Counter()
    bc = whitelist.pop()
    barcode_length = len(bc)
    whitelist.add(bc)

    output_cols = [
        'read_id', 'CR', 'CY', 'UR', 'UY', 'chr', 'start', 'end', 'mapq', 'CB']
    with open(tags_output, 'w') as fh:
        fh.write("\t".join(output_cols))
        fh.write("\n")

    for df_tags in pd.read_csv(barcode_tags, sep='\t', chunksize=chunksize):
        df_tags["CB"] = "-"
        selected = df_tags["CR"].str.len() >= barcode_length - max_ed
        df_tags.loc[selected, "CB"] = df_tags.loc[selected].apply(
            lambda x: determine_barcode(
                x.CR, whitelist, max_ed, min_ed_diff),
            axis=1)
        df_tags[output_cols].to_csv(
            tags_output, mode='a', sep='\t', header=None, index=False)
        barcode_counter.update(df_tags["CB"])

    del barcode_counter["-"]
    return barcode_counter


def main(args):
    """Run main entry point."""
    logger.info("Reading whitelist.")
    whitelist = set(pd.read_csv(
        args.whitelist, index_col=None, sep='\t', header=None)[0])

    logger.info("Processing reads.")
    barcode_counter = process_records(
        args.barcode_tags, whitelist,
        args.max_ed, args.min_ed_diff,
        args.output_tags)

    with open(args.output_counts, "w") as f:
        for bc, n in barcode_counter.most_common():
            f.write(f"{bc}\t{n}\n")
    logger.info("Finished.")
