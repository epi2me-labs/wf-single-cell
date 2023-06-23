#!/usr/bin/env python
"""Extract barcode."""
import collections
from pathlib import Path

import editdistance as ed
import pandas as pd
import parasail
from pysam import AlignmentFile

from .sc_util import kit_adapters  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """
    Parse the command line arguments.

    :return: object containing all supplied arguments
    :rtype: class argparse.Namespace
    """
    # Create argument parser
    parser = wf_parser("extract_barcode")

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing \
            reads aligned to a reference",
        type=Path,
    )

    parser.add_argument(
        "--contig",
        help="Contig/chromosome to process",
    )

    parser.add_argument(
        "superlist",
        help="Comprehensive whitelist of all possible cell barcodes.\
        These vary depending on which 10X kit was used. \
        For 3' v3 single cell gene expression \
        kit: data/3M-february-2018.txt.gz. \
        For 3' v2 single cell gene expression \
        kit: data/737K-august-2016.txt.gz. \
        For 5' single cell gene expression \
        kit: data/737K-august-2016.txt.gz. \
        For single cell multiome (ATAC + GEX) \
        kit: data/737K-arc-v1.txt.gz",
        type=Path,
        default=None,
    )

    # Optional arguments
    parser.add_argument(
        "-k",
        "--kit",
        help="Specify either the 10X 3' gene expression kit (3prime), the 5' \
        gene expression kit (5prime), or the multiome kit (multiome) This \
        determines which adapter sequences to search for in the reads \
        [3prime]",
        default="3prime",
        choices=['3prime', '5prime', 'multiome']
    )

    parser.add_argument(
        "--min_barcode_qv",
        help="Minimum quality score in a barcode for it to be considered \
        a high-quality barcode to be used in whitelist creation [15].",
        default=15,
        type=int,
    )

    parser.add_argument(
        "--adapter1_suff_length",
        help="Use this many suffix bases from adapter1 sequence \
            in the alignment query. For example, specifying 12 \
            would mean that the last 12 bases of the specified \
            read1 sequence will be included in the probe sequence \
            [10]",
        default=10,
        type=int,
    )

    parser.add_argument(
        "-T",
        "--polyt_length",
        help="Length of polyT sequence to use in the \
        alignment query (ignored with --kit=5prime) [10]",
        type=int,
        default=10,
    )

    parser.add_argument(
        "-t",
        "--threads",
        help="Threads to use [4]",
        type=int,
        default=4,
    )

    parser.add_argument(
        "--barcode_length",
        help="Cell barcode length [16]",
        type=int,
        default=16,
    )

    parser.add_argument(
        "--umi_length",
        help="UMI length [12]",
        type=int,
        default=12,
    )

    parser.add_argument(
        "-o", "--gap_open", help="Gap open penalty [2]", type=int, default=2
    )

    parser.add_argument(
        "-e",
        "--gap_extend",
        help="Gap extend penalty [4]",
        type=int,
        default=4)

    parser.add_argument(
        "-m",
        "--match",
        help="Match score [5]",
        type=int,
        default=5)

    parser.add_argument(
        "-x", "--mismatch", help="Mismatch score [-1]", type=int, default=-1
    )

    parser.add_argument(
        "-n",
        "--acg_to_n_match",
        help="Score for A/C/G<-->N match [1]",
        type=int,
        default=1,
    )

    parser.add_argument(
        "-s",
        "--t_to_n_match",
        help="Score for T<-->N match [1]",
        type=int,
        default=1)

    parser.add_argument(
        "-w",
        "--window",
        help="Number of bases to query at start of read [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--max_adapter1_ed",
        help="Max edit distance with the adapter1 sequence (upstream of cell \
        barcode) [3]",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--output_read_tags",
        help="Output TSV file with read_ids and associated tags",
        type=Path
    )

    parser.add_argument(
        "--output_barcode_counts",
        help="Output TSV file containing high-quality barcode counts",
        type=Path
    )

    return parser


def update_matrix(match, mismatch, acg_to_n_match, t_to_n_match):
    """
    Update matrix.

    Create new parasail scoring matrix. 'N' is used as wildcard character
    for barcodes and has its own match parameter (0 per default).

    :param match: Score for matchinging bases
    :param mismatch: Score for mismatches
    :param acg_to_n_match: Score for matching tacq to n
    :param t_to_n_match: Score for matching t to n
    :return: custom parasail alignment matrix
    :rtype: parasail.bindings_v2.Matrix
    """
    matrix = parasail.matrix_create("ACGTN", match, mismatch)

    ############################
    # SCORING MATRIX POSITIONS #
    ############################
    #     A   C   G   T   N
    # A   0   1   2   3   4   5
    # C   6   7   8   9   10  11
    # G   12  13  14  15  16  17
    # T   18  19  20  21  22  23
    # N   24  25  26  27  28  29

    # Update scoring matrix so that N matches A/C/G/N
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = t_to_n_match
    matrix.pointer[0].matrix[27] = t_to_n_match
    return matrix


def ascii_encode_qscores(integers):
    """Convert integer quality values into ASCII characters."""
    return "".join(map(lambda x: chr(x + 33), integers))


def ascii_decode_qscores(string):
    """Convert ASCII character quality values into integers."""
    return list(map(lambda x: ord(x) - 33, string))


def parse_probe_alignment(
        p_alignment, adapter1_probe_seq, barcode_length, umi_length,
        prefix_qual, prefix_seq
        ):
    """Parse probe alignment."""
    ref_alignment = p_alignment.traceback.ref
    query_alignment = p_alignment.traceback.query

    # Find the position of the Ns in the alignment. These correspond
    # to the cell barcode + UMI sequences bound by the read1 and polyT
    bc_start_pos = ref_alignment.find('N')
    if bc_start_pos > -1:
        # umi_end_pos = bc_start_pos + barcode_length + umi_length

        # The alignment of the first N=<args.window> bases of the read
        # to the probe sequence results in the following output, where
        # the alignment around the Ns is anchored by adapter1 and polyT.
        # The positions in the alignment might contain deletions, marked
        # as a dash (-) in the query_alignment:
        #
        #      bc_start_pos               umi_end_pos
        #           v                          v
        # CTTCCGATCTNNNNNNNNNNNNNNNNNNNNNNNNNNNNTTTTTTTTT <-ref_alignment
        # |||||||||||||||||||||||||||||||||||||||||||||||
        # CTTCCGATCT-ATCAGTGATCCACTAGAGGGAGCCGTGTTTTTTTTT <-query_alignment
        # |________||__________________________|
        #  adapter1         barcode_umi

        # The adapter1 sequence comprises the first part of the alignment
        adapter1 = query_alignment[0:bc_start_pos]
        adapter1_ed = ed.eval(adapter1, adapter1_probe_seq)

        # Extract the corresponding positions from both the query sequence
        # and qscores (prefix_qual)
        ascii_q = ascii_encode_qscores(prefix_qual)

        barcode = query_alignment[bc_start_pos: bc_start_pos + barcode_length]
        umi = query_alignment[
            bc_start_pos + barcode_length: bc_start_pos + barcode_length + umi_length]

        barcode_no_ins = barcode.replace("-", "")
        prefix_seq_bc_start = prefix_seq.find(barcode_no_ins)
        prefix_seq_bc_end = prefix_seq_bc_start + len(barcode_no_ins)
        bc_q_ascii = ascii_q[prefix_seq_bc_start:prefix_seq_bc_end]

        umi_no_ins = umi.replace("-", "")
        prefix_seq_umi_start = prefix_seq.find(umi_no_ins)
        prefix_seq_umi_end = prefix_seq_umi_start + len(umi_no_ins)
        umi_q_ascii = ascii_q[prefix_seq_umi_start:prefix_seq_umi_end]

    else:
        # No Ns in the probe successfully aligned -- we will ignore this read
        adapter1_ed = len(adapter1_probe_seq)
        barcode_no_ins = ""
        umi_no_ins = ""
        bc_q_ascii = ""
        umi_q_ascii = ""

    return adapter1_ed, barcode_no_ins, umi_no_ins, bc_q_ascii, umi_q_ascii


def align_adapter(args):
    """Align a single adapter template to read and compute identity.

    :param bam_path: BAM file
    :param chrom: Chromosome/contig to process
    :param: matrix: custom parasail alignment matrix
    :returns: tuple containing two DataFrames
        -  Index:  read_id', Columns: 'CR' 'CY
        -  Index: barcode, Columns: count
    """
    # Build align matrix and define the probe sequence for alignments
    matrix = update_matrix(
        args.match, args.mismatch, args.acg_to_n_match, args.t_to_n_match)

    # Use only the specified suffix length of adapter1
    adapter1_probe_seq = args.adapter1_seq[-args.adapter1_suff_length:]

    if args.kit in ("3prime", "multiome"):
        # Compile the actual probe sequence of
        # <adapter1_suffix>NNN...NNN<TTTTT....>
        probe_seq = "{a1}{bc}{umi}{pt}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            pt="T" * args.polyt_length,
        )
    elif args.kit == "5prime":
        # Compile the actual probe sequence of
        # <adapter1_suffix>NNN...NNN<TTTCTTATATGGG>
        probe_seq = "{a1}{bc}{umi}{tso}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            tso="TTTCTTATATGGG",
        )
    else:
        raise Exception("Invalid kit name! Specify either 3prime or 5prime.")

    with AlignmentFile(
            str(args.bam), "rb") as bam_fh, \
            open(args.output_read_tags, 'w') as tags_fh:

        # Write the header
        tags_fh.write("read_id\tCR\tCY\tUR\tUY\tchr\tstart\tend\tmapq\n")

        barcode_counts = collections.Counter()

        for align in bam_fh.fetch(contig=args.contig):
            if align.is_supplementary:
                continue
            prefix_seq = align.get_forward_sequence()[: args.window]
            prefix_qv = align.get_forward_qualities()[: args.window]

            p_alignment = parasail.sw_trace(
                s1=prefix_seq,
                s2=probe_seq,
                open=args.gap_open,
                extend=args.gap_extend,
                matrix=matrix,
            )

            adapter1_ed, barcode, umi, bc_qscores, umi_qscores \
                = parse_probe_alignment(
                    p_alignment, adapter1_probe_seq, args.barcode_length,
                    args.umi_length, prefix_qv, prefix_seq
                    )

            # Require minimum read1 edit distance
            if (adapter1_ed <= args.max_adapter1_ed) & \
                    (len(barcode) > 0) & (len(umi) > 0):

                bc_min_qv = min(ascii_decode_qscores(bc_qscores))
                if bc_min_qv >= args.min_barcode_qv:
                    barcode_counts[barcode] += 1
                    # nh: I think if we are not including a barcode in the
                    # count due to low min qual, then we should also be
                    # ommiting from the records
                tags_fh.write('\t'.join([
                    align.query_name, barcode, bc_qscores,
                    umi, umi_qscores, args.contig,
                    str(align.get_reference_positions()[0]),
                    str(align.get_reference_positions()[-1]),
                    str(align.mapping_quality)]) + '\n')

    bc_counts = pd.DataFrame.from_dict(
        barcode_counts,
        columns=['count'],
        orient='index').sort_values('count', ascending=False)
    return bc_counts


def main(args):
    """Run entry point."""
    logger = get_named_logger('ExtractBC')
    args.adapter1_seq = kit_adapters[args.kit]['adapter1']
    wl = pd.read_csv(args.superlist, header=None).iloc[:, 0].values

    logger.info(f"Extracting uncorrected barcodes from {args.bam}")

    barcode_counts = align_adapter(args)

    # Filter barcode counts against barcode superlist
    logger.info(
        f"Writing superlist-filtered barcode counts to "
        f"{args.output_barcode_counts}")

    bc_counts = barcode_counts[barcode_counts.index.isin(wl)]
    bc_counts.to_csv(
        args.output_barcode_counts, index=True, sep='\t', header=False)
