"""Extract barcode."""
import collections
from enum import Enum
from pathlib import Path
import sys

import editdistance as ed
import pandas as pd
import parasail
from pysam import FastxFile

from .sc_util import kit_adapters, rev_cmp  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


class KitName(str, Enum):
    """Kit name Enum."""

    prime3 = '3prime'
    prime5 = '5prime'
    multiome = 'multiome'
    visium = 'visium'


def argparser():
    """
    Parse the command line arguments.

    :return: object containing all supplied arguments
    :rtype: class argparse.Namespace
    """
    parser = wf_parser("extract_barcode")

    parser.add_argument(
        "fastq", type=Path,
        help="Read input file.")

    parser.add_argument(
        "superlist", type=Path,
        help="Comprehensive whitelist of all possible cell barcodes.\
        These vary depending on which 10X kit was used.")

    parser.add_argument(
        "--kit", type=KitName, default=KitName.prime3,
        help="Specify the 10x kit name (without version). This \
        determines which adapter sequences to search for in the reads.")

    parser.add_argument(
        "--min_barcode_qv", type=int, default=15,
        help="Minimum quality score in a barcode for it to be considered \
        a high-quality barcode to be used in whitelist creation.")

    parser.add_argument(
        "--adapter1_suff_length", type=int, default=10,
        help="Use this many suffix bases from adapter1 sequence \
            in the alignment query. For example, specifying 12 \
            would mean that the last 12 bases of the specified \
            read1 sequence will be included in the probe sequence.")

    parser.add_argument(
        "--polyt_length", type=int, default=10,
        help="Length of polyT sequence to use in the \
            alignment query (ignored with --kit=5prime).")

    parser.add_argument(
        "--barcode_length", type=int, default=16,
        help="Cell barcode length")

    parser.add_argument(
        "--umi_length", type=int, default=12,
        help="UMI length.")

    parser.add_argument(
        "--gap_open", type=int, default=2,
        help="Gap open penalty.")

    parser.add_argument(
        "--gap_extend", type=int, default=4,
        help="Gap extend penalty.")

    parser.add_argument(
        "--match", type=int, default=5,
        help="Match score.")

    parser.add_argument(
        "--mismatch", type=int, default=-1,
        help="Mismatch score")

    parser.add_argument(
        "--acg_to_n_match", type=int, default=1,
        help="Score for A/C/G<-->N match.")

    parser.add_argument(
        "--t_to_n_match", type=int, default=1,
        help="Score for T<-->N match.")

    parser.add_argument(
        "--window", type=int, default=100,
        help="Number of bases to query at start of read.")

    parser.add_argument(
        "--max_adapter1_ed", type=int, default=3,
        help="Max edit distance with the adapter1 sequence (upstream of cell \
        barcode).")

    # TODO: are these really required arguments?
    parser.add_argument(
        "--output_read_tags",
        help="Output TSV file with read_ids and associated tags"
    )

    parser.add_argument(
        "--output_barcode_counts",
        help="Output TSV file containing high-quality barcode counts"
    )

    return parser


def update_matrix(match, mismatch, acg_to_n_match, t_to_n_match):
    """Update matrix.

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
    """Parse probe alignment.

    Regardless of kit, the sequences need to be in adapter orientation.
    """
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

        barcode = query_alignment[bc_start_pos: bc_start_pos + barcode_length]
        umi = query_alignment[
            bc_start_pos + barcode_length: bc_start_pos + barcode_length + umi_length]

        barcode_no_ins = barcode.replace("-", "")
        prefix_seq_bc_start = prefix_seq.find(barcode_no_ins)
        prefix_seq_bc_end = prefix_seq_bc_start + len(barcode_no_ins)
        bc_q_ascii = prefix_qual[prefix_seq_bc_start:prefix_seq_bc_end]

        umi_no_ins = umi.replace("-", "")
        prefix_seq_umi_start = prefix_seq.find(umi_no_ins)
        prefix_seq_umi_end = prefix_seq_umi_start + len(umi_no_ins)
        umi_q_ascii = prefix_qual[prefix_seq_umi_start:prefix_seq_umi_end]

    else:
        # No Ns in the probe successfully aligned -- we will ignore this read
        adapter1_ed = len(adapter1_probe_seq)
        barcode_no_ins = ""
        umi_no_ins = ""
        bc_q_ascii = ""
        umi_q_ascii = ""

    return (
        adapter1_ed, barcode_no_ins, umi_no_ins,
        bc_q_ascii, umi_q_ascii
    )


def align_adapter(args, fastq_out=sys.stdout):
    """
    Align a single adapter template to read and compute identity.

    :param fastq: FASTQ file
    :param chrom: Chromosome/contig to process
    :param: matrix: custom parasail alignment matrix
    :returns: tuple containing two DataFrames
        -  Index:  read_id', Columns: 'CR' 'CY
        -  Index: barcode, Columns: count
    """
    # Build align matrix and define the probe sequence for alignments
    # Note: once CW-2853 is done, the 3prime reads will need flipping
    matrix = update_matrix(
        args.match, args.mismatch, args.acg_to_n_match, args.t_to_n_match)

    # Use only the specified suffix length of adapter1
    adapter1_probe_seq = args.adapter1_seq[-args.adapter1_suff_length:]

    if args.kit in (KitName.prime3, KitName.multiome, KitName.visium):
        # For these kits the probe needs to be reverse complemented
        # <adapter1_suffix>NNN...NNN<TTTTT....>
        probe_seq = "{a1}{bc}{umi}{pt}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            pt="T" * args.polyt_length
        )
    elif args.kit == KitName.prime5:
        # Compile the actual probe sequence of
        # <adapter1_suffix>NNN...NNN<TTTCTTATATGGG>
        probe_seq = "{a1}{bc}{umi}{tso}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            tso="TTTCTTATATGGG"
        )
    else:
        raise Exception(
            "Invalid kit_name parameter! Specify either 3prime, multiome, or 5prime.")

    with FastxFile(
            str(args.fastq), "rb") as fastq_fh, \
            open(args.output_read_tags, 'w') as tags_fh:

        # Write the header
        tags_fh.write("read_id\tCR\tCY\tUR\tUY\n")

        barcode_counts = collections.Counter()

        for read in fastq_fh:

            if args.kit in (KitName.prime3, KitName.multiome, KitName.visium):
                # Flip back to barcode orientation (reverse)
                prefix_seq = rev_cmp(read.sequence)[: args.window]
                prefix_qv = read.quality[::-1][: args.window:]
            else:
                prefix_seq = read.sequence[: args.window]
                prefix_qv = read.quality[: args.window]

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

                # Escape double quotes, with a preceding `"` in the quality strings
                # see https://rfc-editor.org/rfc/rfc4180.html
                umi_q_quoted = umi_qscores.replace('"', '""')
                barcode_q_quoted = bc_qscores.replace('"', '""')

                tags_fh.write('\t'.join([
                    f'"{read.name}"', f'"{barcode}"', f'"{barcode_q_quoted}"',
                    f'"{umi}"', f'"{umi_q_quoted}"']) + '\n')
                # For full length reads, adapter2 already trimmed.
                # Now the barcode and UMI have been extracted, these along with adapter1
                # can be removed.
                if args.kit in (KitName.prime3, KitName.multiome, KitName.visium):
                    # The reads will be cDNA-polyA-UMI-BC-Adapter1,
                    # so to trim from right
                    trim_side = 'right'
                    trim_pos = p_alignment.end_query - 1
                    trim_pos -= args.polyt_length
                else:
                    # 5prime kit
                    # The reads will be adatper1-BC-UMI-TSO-cDNA-polyT,
                    # so trim from left
                    trim_side = 'left'
                    trim_pos = p_alignment.end_query

                # trim and write the read
                if trim_side == 'right':
                    read.sequence = read.sequence[:-trim_pos]
                    read.quality = read.quality[:-trim_pos]
                else:
                    read.sequence = read.sequence[trim_pos:]
                    read.quality = read.quality[trim_pos:]
                if read.sequence:
                    fastq_out.write(str(read) + '\n')

    bc_counts = pd.DataFrame.from_dict(
        barcode_counts,
        columns=['count'],
        orient='index').sort_values('count', ascending=False)
    bc_counts.index.name = 'barcode'
    return bc_counts


def main(args):
    """Run entry point."""
    logger = get_named_logger('ExtractBC')
    args.adapter1_seq = kit_adapters[args.kit]['adapter1']

    logger.info(f"Loading barcode whitelist from {args.superlist}")
    wl = pd.read_csv(args.superlist, header=None).iloc[:, 0].values

    logger.info(f"Extracting uncorrected barcodes from {args.fastq}")
    barcode_counts = align_adapter(args)

    logger.info(
        f"Writing superlist-filtered barcode counts to "
        f"{args.output_barcode_counts}")
    bc_counts = barcode_counts[barcode_counts.index.isin(wl)]
    bc_counts.to_csv(
        args.output_barcode_counts, index=True, sep='\t')

    logger.info("Finished.")
