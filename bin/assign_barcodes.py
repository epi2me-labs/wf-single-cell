import argparse
import collections
import gzip
import logging
import math
import multiprocessing
import os
import pathlib
import re
import shutil
import tempfile

import editdistance as ed
import parasail
import pysam
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing reads aligned to a reference. \
            Alignments must have the CR and CY tags.",
        type=str,
    )

    parser.add_argument(
        "whitelist", help="File containing list of expected cell barcodes", type=str
    )

    # Optional arguments
    parser.add_argument(
        "-t", "--threads", help="Threads to use [4]", type=int, default=4
    )

    parser.add_argument(
        "-k", help="Kmer size to use for whitelist filtering [5]", type=int, default=5
    )

    parser.add_argument(
        "--output_bam",
        help="Output BAM file containing aligned reads with tags for uncorrected \
        barcodes (CR), corrected barcodes (CB), barcode QVs (CY), uncorrected \
        UMIs (UR), and UMI QVs (UY) [bc_corr.umi_uncorr.sorted.bam]",
        type=str,
        default="bc_corr.umi_uncorr.sorted.bam",
    )

    parser.add_argument(
        "--output_counts",
        help="Output TSV file containing counts for each of the assigned \
        barcodes [barcode_counts.tsv]",
        type=str,
        default="barcode_counts.tsv",
    )

    parser.add_argument(
        "--max_ed",
        help="Max edit distance between putative barcode \
                        and the matching whitelist barcode [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--min_ed_diff",
        help="Min difference in edit distance between the \
                        (1) putative barcode vs top hit and (2) putative \
                        barcode vs runner-up hit [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--kit",
        help="Specify either the 10X 3' gene expression kit (3prime), the 5' \
        gene expression kit (5prime), or the multiome kit (multiome) This \
        determines which adapter sequences to search for in the reads \
        [3prime]",
        default="3prime",
        type=str,
    )

    parser.add_argument(
        "--adapter1_suff_length",
        help="Use this many suffix bases from adapter1 sequence  in the \
            alignment query. For example, specifying 12 would mean that the last \
            12 bases of the specified read1 sequence will be included in the \
            probe sequence [10]",
        default=10,
        type=int,
    )

    parser.add_argument(
        "-T",
        "--polyT_length",
        help="Length of polyT sequence to use in the alignment query (ignored \
        with --kit=5prime) [10]",
        type=int,
        default=10,
    )

    parser.add_argument(
        "--barcode_length", help="Cell barcode length [16]", type=int, default=16
    )

    parser.add_argument("--umi_length", help="UMI length [12]", type=int, default=12)

    parser.add_argument(
        "-w",
        "--window",
        help="Number of bases to query at start of read [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "-o", "--gap_open", help="Gap open penalty [2]", type=int, default=2
    )

    parser.add_argument(
        "-e", "--gap_extend", help="Gap extend penalty [4]", type=int, default=4
    )

    parser.add_argument("-m", "--match", help="Match score [5]", type=int, default=5)

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
        "-s", "--t_to_n_match", help="Score for T<-->N match [1]", type=int, default=1
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    # verify kit selection
    if (args.kit != "3prime") and (args.kit != "5prime") and (args.kit != "multiome"):
        raise Exception(
            "Invalid kit name! Specify either 3prime, 5prime or \
        multiome."
        )

    if (args.kit == "3prime") or (args.kit == "multiome"):
        # Read1 adapter
        args.adapter1_seq = "CTACACGACGCTCTTCCGATCT"
        # TSO adapter
        args.adapter2_seq = "ATGTACTCTGCGTTGATACCACTGCTT"
    elif args.kit == "5prime":
        # Read1 adapter
        args.adapter1_seq = "CTACACGACGCTCTTCCGATCT"
        # Poly-dT RT adapter
        args.adapter2_seq = "GTACTCTGCGTTGATACCACTGCTT"

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_bam)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def find(target, my_string):
    """
    Return indices corresponding to the positions of <target> in <my_string>.

    :param target: character to search for in <my_string>
    :type target: str
    :param my_string: string to search for <target>
    :type my_string: str
    :return: generator returning indices corresponding to the position(s) of
        <target> in <my_string>
    :rtype: int
    """
    for i in range(len(my_string)):
        if my_string[i] == target:
            yield i


def update_matrix(args):
    """
    Create new parasail scoring matrix. 'N' is used as wildcard character
    for barcodes and has its own match parameter (0 per default).

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: custom parasail alignment matrix
    :rtype: parasail.bindings_v2.Matrix
    """
    matrix = parasail.matrix_create("ACGTN", args.match, args.mismatch)

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
    # pointers = [4, 9, 14, 20, 21, 22, 24]
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = args.acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = args.t_to_n_match
    matrix.pointer[0].matrix[27] = args.t_to_n_match
    return matrix


def calc_ed_with_whitelist(bc_uncorr, whitelist):
    """
    Find minimum and runner-up barcode edit distance by iterating through the
    whitelist of expected barcodes.

    :param bc_uncorr: Uncorrected cell barcode
    :type bc_uncorr: str
    :param whitelist: Filtered whitelist of cell barcodes
    :type whitelist: list
    :return: Corrected barcode assignment, edit distance, and difference in edit
        distance between the top match and the next closest match
    :rtype: str, int, int
    """
    bc_match = "X" * len(bc_uncorr)
    bc_match_ed = len(bc_uncorr)
    next_bc_match_ed = len(bc_uncorr)
    for wl_bc in whitelist:
        d = ed.eval(bc_uncorr, wl_bc)
        if d < bc_match_ed:
            next_bc_match_ed = bc_match_ed
            bc_match_ed = d
            bc_match = wl_bc
        elif d < next_bc_match_ed:
            next_bc_match_ed = d
    next_match_diff = next_bc_match_ed - bc_match_ed

    return bc_match, bc_match_ed, next_match_diff


LOOKUP = []

for q in range(100):
    LOOKUP.append(pow(10, -0.1 * q))


def compute_mean_qscore(scores):
    """
    Returns the phred score corresponding to the mean of the probabilities
    associated with the phred scores provided.

    :param scores: Iterable of phred scores.
    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if not scores:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += LOOKUP[val]
    mean_prob = sum_prob / len(scores)

    return -10.0 * math.log10(mean_prob)


def find_feature_qscores(feature, p_alignment, prefix_seq, prefix_qv):
    """
    Using the parasail alignment results, find the qscores corresponding to the
    feature (e.g. barcode or UMI) positions in the read.

    :param feature: Feature sequence identified from the parasail alignment
    :type feature: str
    :param p_alignment: Parasail alignment object
    :type p_alignment: class 'parasail.bindings_v2.Result'
    :param prefix_seq: Nucleotide sequence from the first <args.window> bp of
        the read
    :type prefix_seq: str
    :param prefix_qv: Qscores from the first <args.window> bp of the read
    :type prefix_qv: np.array
    :return: Array of phred scale qscores from the identified feature region
    :rtype: np.array
    """
    # Strip feature alignment string of insertions (-)
    feature_no_ins = feature.replace("-", "")

    # Find where the stripped feature starts and ends in the larger prefis_seq
    prefix_seq_feature_start = prefix_seq.find(feature_no_ins)
    prefix_seq_feature_end = prefix_seq_feature_start + len(feature_no_ins)

    # Use these start/end indices to locate the correspoding qscores in prefix_qv
    feature_qv = prefix_qv[prefix_seq_feature_start:prefix_seq_feature_end]
    feature_qv_ascii = "".join(map(lambda x: chr(x + 33), feature_qv))

    return feature_qv_ascii


def parse_probe_alignment(p_alignment, align, prefix_seq, prefix_qv):
    """
    Parse a parasail alignment alignment and add uncorrected UMI and UMI QV
    values as tags to the BAM alignment.

    :param p_alignment: parasail alignment
    :type p_alignment: class 'parasail.bindings_v2.Result'
    :param align: pysam BAM alignment
    :type align: class 'pysam.libcalignedsegment.AlignedSegment'
    :return: pysam BAM alignment with the UR and UY tags added
    :rtype: class 'pysam.libcalignedsegment.AlignedSegment'
    """
    # Find the position of the Ns in the parasail alignment. These correspond
    # to the UMI sequences bound by the cell barcode and polyT
    idxs = list(find("N", p_alignment.traceback.ref))
    if len(idxs) > 0:
        umi = p_alignment.traceback.query[min(idxs) : max(idxs) + 1]

        qscores = find_feature_qscores(umi, p_alignment, prefix_seq, prefix_qv)

        # print(p_alignment.traceback.ref)
        # print(p_alignment.traceback.comp)
        # print(p_alignment.traceback.query)
        # print()

        umi = umi.replace("-", "")

        # Uncorrected UMI = UR:Z
        align.set_tag("UR", umi, value_type="Z")
        # UMI quality score = UY:Z
        align.set_tag("UY", qscores, value_type="Z")

    return align


def get_uncorrected_umi(align, args):
    """
    Aligns a probe sequence containing the adapter1+corrected_barcode+Ns+polyT to
    the read. Bases aligning to the Ns in the probe sequence correspond to the
    UMI positions. Extract those bases and consider those to be the uncorrected
    UMI sequence.

    :param align: pysam BAM alignment with the CB tag
    :type align: class 'pysam.libcalignedsegment.AlignedSegment'
    :param args: object containing all supplied arguments
    :type args: class 'argparse.Namespace'
    :return: pysam BAM alignment with the UR and UY tags added
    :rtype: class 'pysam.libcalignedsegment.AlignedSegment'
    """
    prefix_seq = align.get_forward_sequence()[: args.window]
    prefix_qv = align.get_forward_qualities()[: args.window]

    # Use only the specified suffix length of adapter1
    adapter1_probe_seq = args.adapter1_seq[-args.adapter1_suff_length :]

    if (args.kit == "3prime") or (args.kit == "multiome"):
        # Compile the actual query sequence of <adapter1_suffix><bc_corr>NNN...N<TTTTT....>
        probe_seq = "{r}{bc}{umi}{pT}".format(
            r=adapter1_probe_seq,
            bc=align.get_tag("CB"),
            umi="N" * args.umi_length,
            pT="T" * args.polyT_length,
        )
    elif args.kit == "5prime":
        # Compile the actual probe sequence of <adapter1_suffix>NNN...NNN<TTTCTTATATGGG>
        probe_seq = "{a1}{bc}{umi}{tso}".format(
            a1=adapter1_probe_seq,
            bc=align.get_tag("CB"),
            umi="N" * args.umi_length,
            tso="TTTCTTATATGGG",
        )
    else:
        raise Exception("Invalid kit name! Specify either 3prime or 5prime.")

    matrix = update_matrix(args)
    parasail_alg = parasail.sw_trace

    p_alignment = parasail_alg(
        s1=prefix_seq,
        s2=probe_seq,
        open=args.gap_open,
        extend=args.gap_extend,
        matrix=matrix,
    )

    # print(p_alignment.traceback.ref)
    # print(p_alignment.traceback.comp)
    # print(p_alignment.traceback.query)
    # print()

    align = parse_probe_alignment(p_alignment, align, prefix_seq, prefix_qv)

    return align


def process_bam_records(tup):
    """
    Process BAM records to assign each read a corrected cell barcode and an
    uncorrected UMI. Do this by loading and processing the barcode whitelist
    then iterating over alignments.
    For each alignment:
    1. Calculate edit distance between uncorrected barcode and barcodes in whitelist
    2.

    :param tup: Tuple containing the input arguments
    :type tup: tup
    :return: Path to a temporary BAM file
    :rtype: str
    """
    input_bam = tup[0]
    chrom = tup[1]
    args = tup[2]

    # Load barcode whitelist and map kmers to indices in whitelist for faster
    # barcode matching
    whitelist, kmer_to_bc_index = load_whitelist(args.whitelist, args.k)

    # Open input BAM file
    bam = pysam.AlignmentFile(input_bam, "rb")

    # Write temp file or straight to output file depending on use case
    if args.threads > 1:
        # Open temporary output BAM file for writing
        suff = f".{chrom}.bam"
        chrom_bam = tempfile.NamedTemporaryFile(
            prefix="tmp.align.", suffix=suff, dir=args.tempdir, delete=False
        )
        bam_out_fn = chrom_bam.name
    else:
        bam_out_fn = args.output_bam

    bam_out = pysam.AlignmentFile(args.output_bam, "wb", template=bam)

    barcode_counter = collections.Counter()

    for align in bam.fetch(contig=chrom):
        # Make sure each alignment in this BAM has an uncorrected barcode and
        # barcode QV
        assert align.has_tag("CR") and align.has_tag("CY"), "CR or CY tags not found"

        bc_uncorr = align.get_tag("CR")

        # Don't consider any uncorrected barcodes that are shorter than k
        if len(bc_uncorr) >= args.k:
            # Decompose uncorrected barcode into N k-mers
            bc_uncorr_kmers = split_seq_into_kmers(bc_uncorr, args.k)
            # Filter the whitelist to only those with at least one of the k-mers
            # from the uncorrected barcode
            filt_whitelist = filter_whitelist_by_kmers(
                whitelist, bc_uncorr_kmers, kmer_to_bc_index
            )

            # Calc edit distances between uncorrected barcode and the filtered
            # whitelist barcodes
            bc_match, bc_match_ed, next_match_diff = calc_ed_with_whitelist(
                bc_uncorr, filt_whitelist
            )

            # Check barcode match edit distance and difference to runner-up edit distance
            condition1 = bc_match_ed <= args.max_ed
            condition2 = next_match_diff >= args.min_ed_diff
            if condition1 and condition2:
                # Add corrected cell barcode = CB:Z
                align.set_tag("CB", bc_match, value_type="Z")

                # Add corrected barcode to probe sequence to fish out uncorrected UMI
                align = get_uncorrected_umi(align, args)

            # Only write BAM entry in output file if we've assigned a corrected
            # barcode and an uncorrected UMI
            if align.has_tag("CB") and align.has_tag("UR"):
                bam_out.write(align)
                barcode_counter[bc_match] += 1

    bam.close()
    bam_out.close()

    return bam_out_fn, barcode_counter


def launch_pool(func, func_args, procs=1):
    """
    Use multiprocessing library to create pool and map function calls to
    that pool

    :param procs: Number of processes to use for pool
    :type procs: int, optional
    :param func: Function to exececute in the pool
    :type func: function
    :param func_args: List containing arguments for each call to function <funct>
    :type func_args: list
    :return: List of results returned by each call to function <funct>
    :rtype: list
    """
    p = multiprocessing.Pool(processes=procs)
    try:
        results = list(tqdm(p.imap(func, func_args), total=len(func_args)))
        p.close()
        p.join()
    except KeyboardInterrupt:
        p.terminate()
    return results


def filter_whitelist_by_kmers(wl, kmers, kmer_to_bc_index):
    """
    Given a list of whitelisted barcodes, return just the
    subset that contain any of the kmers contained in the
    query barcode.

    :param wl: Full barcode whitelist
    :type wl: list
    :param kmers: K-mers to use for whitelist filtering
    :type kmers: list
    :param kmer_to_bc_index: Map of k-mers to the whitelist indices corresponding
        to all barcodes containing that k-mer
    :type kmer_to_bc_index: dict
    :return: List of filtered barcodes
    :rtype: list
    """
    # collect sets of indices that each kmer points to
    id_sets = [
        kmer_to_bc_index[kmer] for kmer in kmers if kmer in kmer_to_bc_index.keys()
    ]

    # retain all barcodes that have at least one kmer match with the query barcode
    all_filt_indices = list(set().union(*id_sets))
    filt_wl = [wl[i] for i in all_filt_indices]
    return filt_wl


def split_seq_into_kmers(seq, k):
    """
    Decompose the supplied <seq> into N=len(seq)-k+1 k-mers.

    :param seq: String of nucleotides
    :type seq: str
    :param k: k-mer length
    :type k: int
    :return: List of k-mers
    :rtype: list
    """
    assert len(seq) >= k, "Pick a value for k that is less than len(barcode)"

    kmers = []
    for i in range(0, len(seq) - k + 1):
        kmer = seq[i : i + k]
        kmers.append(kmer)
    return kmers


def load_whitelist(whitelist, k=5):
    """
    Read in barcode whitelist and create dictionary mapping each k-mer to all
    barcodes in the whitelist containing that k-mer.

    :param whitelist: Path to the barcode whitelist
    :type whitelist: str
    :param k: k-mer length
    :type k: int
    :return: List of whitelisted barcodes and dictionary mapping all k-mers to
        indices in the whitelist corresponding to barcodes containing that k-mer
    :rtype: list, dict
    """
    wl = []
    with open(whitelist) as file:
        for line in file:
            bc = line.strip().split("-")[0]
            wl.append(bc)

    wl.sort()
    kmer_to_bc_index = {}
    for index, bc in enumerate(wl):
        bc_kmers = split_seq_into_kmers(bc, k)
        for bc_kmer in bc_kmers:
            if bc_kmer not in kmer_to_bc_index.keys():
                kmer_to_bc_index[bc_kmer] = set([index])
            else:
                kmer_to_bc_index[bc_kmer].add(index)
    return wl, kmer_to_bc_index


def get_bam_info(bam):
    """
    Use `samtools idxstat` to get number of alignments and names of all contigs
    in the reference.

    :param bam: Path to sorted BAM file
    :type bame: str
    :return: Sum of all alignments in the BAM index file and list of all chroms
    :rtype: int,list
    """
    bam = pysam.AlignmentFile(bam, "rb")
    stats = bam.get_index_statistics()
    n_aligns = int(sum([contig.mapped for contig in stats]))
    chroms = dict(
        [(contig.contig, contig.mapped) for contig in stats if contig.mapped > 0]
    )
    bam.close()
    return n_aligns, chroms


def main(args):
    init_logger(args)
    # logger.info("Getting BAM statistics")
    n_reads, chroms = get_bam_info(args.bam)

    if args.threads > 1:
        # Create temporary directory
        if os.path.exists(args.tempdir):
            shutil.rmtree(args.tempdir, ignore_errors=True)
        os.mkdir(args.tempdir)

        # Process BAM alignments from each chrom separately
        logger.info(f"Assigning barcodes to reads in {args.bam}")
        func_args = []
        chroms_sorted = dict(sorted(chroms.items(), key=lambda item: item[1]))
        for chrom in chroms_sorted.keys():
            func_args.append((args.bam, chrom, args))

        results = launch_pool(process_bam_records, func_args, args.threads)
        chrom_bam_fns, barcode_counters = zip(*results)

        barcode_counter = sum(barcode_counters, collections.Counter())

        tmp_bam = tempfile.NamedTemporaryFile(
            prefix="tmp.align.", suffix=".unsorted.bam", dir=args.tempdir, delete=False
        )
        merge_parameters = ["-f", tmp_bam.name] + list(chrom_bam_fns)
        pysam.merge(*merge_parameters)

        pysam.sort("-@", str(args.threads), "-o", args.output_bam, tmp_bam.name)

        logger.info("Cleaning up temporary files")
        shutil.rmtree(args.tempdir, ignore_errors=True)

    else:
        # Hopefully the chromosome name is prefix of BAM filename
        REGEX = r"([A-Za-z0-9.]+).sorted.bam"
        m = re.search(REGEX, args.bam)
        chrom = m.group(1)
        func_args = (args.bam, chrom, args)
        chrom_bam_fn, barcode_counter = process_bam_records(func_args)

    with open(args.output_counts, "w") as f:
        for bc, n in barcode_counter.most_common():
            f.write(f"{bc}\t{n}\n")


if __name__ == "__main__":
    args = parse_args()

    main(args)
