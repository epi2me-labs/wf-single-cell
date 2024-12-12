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
        "report", type=Path,
        help="Path to TSV file to store reasons for barcode assignment.")

    parser.add_argument(
        "--chunksize", type=int, default=50000,
        help="Process the BAM in chunks no larger than this.")

    parser.add_argument(
        "--use_kmer_index", action='store_true',
        help="Use a kmer index to reduce the search space of fuzzy matching.")

    parser.add_argument(
        "--max_ed", type=int, default=2,
        help="Max. edit distance between putative barcode \
            and the matching whitelist barcode.")

    parser.add_argument(
        "--min_ed_diff", type=int, default=2,
        help="Min. difference in edit distance between the \
            best and second best whitelist matches.")

    return parser


def determine_barcode(
        bc_uncorr, whitelist, whiteset,
        max_ed, min_ed_diff, assignment_log, index=None):
    """Find barcode in a whitelist corresponding to read barcode.

    :param bc_uncorr: uncorrected barcode.
    :param whitelist: list of possible barcodes.
    :param whiteset: whitelist as a set.
    :param max_ed: max. edit distance between barcode and whitelist hit.
    :param min_ed_diff: min. edit distance difference between first and
        second best hits in order to accept the first as valid.
    :param assignment_log: a Counter object to store reasons for barcode assignment.
    :param index: a kmer index for reducing search space of fuzzy-matching.

    Passing the whitelist as both a list and set is for performance reasons
    when calling this function many times.
    """
    # quick return
    if bc_uncorr in whiteset:
        assignment_log["bc_shortlist_exact_match"] += 1
        return bc_uncorr

    if index is not None:
        shortlist = set()
        for kmer in build_kmers(bc_uncorr):
            shortlist.update(index[kmer])
        shortlist = list(shortlist)
    else:
        shortlist = whitelist

    result = extract(
        bc_uncorr,
        shortlist,
        scorer=rapidfuzz.distance.Levenshtein.distance,
        score_cutoff=max_ed + min_ed_diff + 1)

    corrected = "-"
    if len(result) > 0:
        bc_match = result[0][0]
        bc_match_ed = result[0][1]
    else:
        assignment_log['bc_no_shortlist_match'] += 1
        return corrected
    if len(result) > 1:
        next_match_diff = result[1][1] - bc_match_ed
    else:
        next_match_diff = len(bc_uncorr)

    # are we better than the second place?
    # This criteria is a little odd: we have (2, 2) as the defaults
    # for max_ed and min_ed_diff. But some true barcodes are within
    # and edit distance of 2 to start with, so they would be guaranteed
    # to be filtered out (the exact match shortcut above saves us a lot
    # of the time). Consider removing this?
    if (bc_match_ed <= max_ed) and (next_match_diff >= min_ed_diff):
        corrected = bc_match
        assignment_log['bc_corrected'] += 1
    else:
        assignment_log['bc_shortlist_multiple_hits'] += 1

    return corrected


def build_index(whitelist, klen=5):
    """Build a kmer index of a list of sequences."""
    index = collections.defaultdict(set)
    for seq in whitelist:
        for ss in build_kmers(seq, klen):
            index[ss].add(seq)
    return index


def build_kmers(seq, klen=5):
    """Create a list of kmers in a sequence."""
    return [seq[i:i+klen] for i in range(0, len(seq) - klen)]


def process_records(
        barcode_tags, whiteset, max_ed, min_ed_diff, tags_output,
        chunksize=50000, use_kmer_index=False):
    """Process read barcodes stored in text file to find whitelist equivalents.

    :param barcode_tags: path to TSV with tag data
    :param whiteset: set of allowed barcodes.
    :param: max_ed: max allowed edit distance between an uncorrected barcode
        and a potential corrected whiteset barcode.
    :param: min_ed_diff: minimum allowed edit distance between top two
        barcode candidates.
    """
    barcode_counter = collections.Counter()
    # we need a list for indexing and because rapidfuzz appears to coerce
    # its input to a list on every call, saves 10% of the time.
    whitelist = list(whiteset)
    barcode_length = len(whitelist[0])
    # for 16mers with 2 mismatches we must have a least a 5mer match.
    # The limit is reached by distributing the mismatches evenly, any
    # perturbation will increase the longest match length.
    # 0123456789ABCDEF
    #     |     |
    index = None
    if use_kmer_index:
        kmer = barcode_length // (max_ed + 1)
        index = build_index(whitelist, klen=kmer)

    output_cols = [
        'read_id', 'CR', 'CY', 'UR', 'UY', 'chr', 'start', 'end', 'mapq', 'CB']
    with open(tags_output, 'w') as fh:
        fh.write("\t".join(output_cols))
        fh.write("\n")

    total_reads = 0
    assignment_log = collections.Counter()
    for df_tags in pd.read_csv(barcode_tags, sep='\t', chunksize=chunksize):
        df_tags["CB"] = "-"
        selected = df_tags["CR"].str.len() >= barcode_length - max_ed
        df_tags.loc[selected, "CB"] = df_tags.loc[selected].apply(
            lambda x: determine_barcode(
                x.CR, whitelist, whiteset, max_ed, min_ed_diff, assignment_log, index),
            axis=1)
        df_tags[output_cols].to_csv(
            tags_output, mode='a', sep='\t', header=None, index=False)
        barcode_counter.update(df_tags["CB"])
        total_reads += len(df_tags)
        logger.info(f"Processed {total_reads} reads.")

    del barcode_counter["-"]
    return barcode_counter, assignment_log


def main(args):
    """Run main entry point."""
    logger.info("Reading whitelist.")
    whiteset = set(pd.read_csv(
        args.whitelist, index_col=None, sep='\t', header=None)[0])

    logger.info("Processing reads.")
    barcode_counter, assignment_log = process_records(
        args.barcode_tags, whiteset,
        args.max_ed, args.min_ed_diff,
        args.output_tags,
        chunksize=args.chunksize,
        use_kmer_index=args.use_kmer_index)

    with open(args.output_counts, "w") as f:
        for bc, n in barcode_counter.most_common():
            f.write(f"{bc}\t{n}\n")

    df_summary = pd.DataFrame.from_dict(assignment_log, orient='index')
    df_summary.to_csv(args.report, sep='\t', header=False)

    logger.info("Finished.")
