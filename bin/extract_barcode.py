import argparse
import collections
import gzip
import logging
import math
import mmap
import multiprocessing
import os
import pathlib
import shutil
import sys
import tempfile

import editdistance as ed
import parasail
import pysam
from tqdm import tqdm

# Make logger globally accessible to all functions
logger = logging.getLogger(__name__)


def parse_args():
    """
    Parse the command line arguments

    :return: object containing all supplied arguments
    :rtype: class argparse.Namespace
    """
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bam",
        help="Sorted BAM file of stranded sequencing reads aligned to a reference",
        type=str,
    )

    parser.add_argument(
        "superlist",
        help="Comprehensive whitelist of all possible cell barcodes. These vary \
        depending on which 10X kit was used. For 3' single cell gene expression \
        kit: data/3M-february-2018.txt.gz. For 5' single cell gene expression \
        kit: data/737K-august-2016.txt. For single cell multiome (ATAC + GEX) \
        kit: data/737K-arc-v1.txt.gz",
        type=str,
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
        type=str,
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
        "--polyT_length",
        help="Length of polyT sequence to use in the alignment query (ignored \
        with --kit=5prime) [10]",
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
        "--output_bam",
        help="Output BAM file containing aligned reads with tags for uncorrected \
        barcodes (CR) and barcode QVs (CY) [bc_uncorr.sorted.bam]",
        type=str,
        default="bc_uncorr.sorted.bam",
    )

    parser.add_argument(
        "--output_barcodes",
        help="Output TSV file containing high-quality barcode counts \
        [barcodes_counts.tsv]",
        type=str,
        default="barcodes_counts.tsv",
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
    elif args.kit == "5prime":
        # Read1 adapter
        args.adapter1_seq = "CTACACGACGCTCTTCCGATCT"

    # Create temp dir and add that to the args object
    p = pathlib.Path(args.output_bam)
    output_dir = p.parents[0]
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=output_dir)
    args.tempdir = tempdir.name

    return args


def init_logger(args):
    """
    Initialize the logger using the specified verbosity level.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


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
    pointers = [4, 10, 16, 24, 25, 26]
    for i in pointers:
        matrix.pointer[0].matrix[i] = args.acg_to_n_match

    # Update N to T matching score so that we start
    # aligning dT sequence instead of Ns.
    matrix.pointer[0].matrix[22] = args.t_to_n_match
    matrix.pointer[0].matrix[27] = args.t_to_n_match
    return matrix


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


def find(char, string):
    """
    Return iterator of indices for positions in a string
    corresponding to a target character

    :param char: Target character whose positions you want to locate in the string
    :type char: str
    :param string: String to search for the target character positions
    :type string: str
    :return: Indices in the string corresponding to the target character
    :rtype: iterator
    """
    for i in range(len(string)):
        if string[i] == char:
            yield i


def edit_distance(query, target):
    """
    Return Levenshtein distance between the two supplied strings

    :param query: Query string to compare against the target
    :type query: str
    :param target: Target string to compare against the query
    :type target: str
    :return: Calculated Levenshtein distance between query and target
    :rtype: int
    """
    d = ed.eval(query, target)
    return d


def parse_probe_alignment(p_alignment, adapter1_probe_seq, args):
    """ """
    # Find the position of the Ns in the alignment. These correspond
    # to the cell barcode + UMI sequences bound by the read1 and polyT
    idxs = list(find("N", p_alignment.traceback.ref))
    if len(idxs) > 0:
        # The Ns in the probe successfully aligned to sequence
        bc_start = min(idxs)

        # The read1 adapter comprises the first part of the alignment
        adapter1 = p_alignment.traceback.query[0:bc_start]
        adapter1_ed = edit_distance(adapter1, adapter1_probe_seq)

        # The barcode + UMI sequences in the read correspond to the
        # positions of the aligned Ns in the probe sequence
        barcode = p_alignment.traceback.query[
            bc_start : (bc_start + args.barcode_length)
        ]
    else:
        # No Ns in the probe successfully aligned -- we will ignore this read
        adapter1_ed = len(adapter1_probe_seq)
        barcode = ""
        bc_start = 0

    return adapter1_ed, barcode, bc_start


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


def align_adapter(tup):
    """
    Aligns a single adapter template to the read an computes the
    identity for the alignment

    :param tup: Tuple containing the function arguments
    :type tup: tup
    :return: Path to temporary BAM containing CR and CY tags, plus a counter
        tracking the number of each barcode that we encounter
    :rtype: str, class 'collections.Counter'
    """

    bam_path = tup[0]
    chrom = tup[1]
    args = tup[2]

    # Build align matrix and define the probe sequence for alignments
    matrix = update_matrix(args)
    parasail_alg = parasail.sw_trace

    # Use only the specified suffix length of adapter1
    adapter1_probe_seq = args.adapter1_seq[-args.adapter1_suff_length :]

    if (args.kit == "3prime") or (args.kit == "multiome"):
        # Compile the actual probe sequence of <adapter1_suffix>NNN...NNN<TTTTT....>
        probe_seq = "{a1}{bc}{umi}{pT}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            pT="T" * args.polyT_length,
        )
    elif args.kit == "5prime":
        # Compile the actual probe sequence of <adapter1_suffix>NNN...NNN<TTTCTTATATGGG>
        probe_seq = "{a1}{bc}{umi}{tso}".format(
            a1=adapter1_probe_seq,
            bc="N" * args.barcode_length,
            umi="N" * args.umi_length,
            tso="TTTCTTATATGGG",
        )
    else:
        raise Exception("Invalid kit name! Specify either 3prime or 5prime.")

    bam = pysam.AlignmentFile(bam_path, "rb")

    # Write output BAM file
    suff = f".{chrom}.bam"
    chrom_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=suff, dir=args.tempdir, delete=False
    )
    bam_out = pysam.AlignmentFile(chrom_bam.name, "wb", template=bam)

    chrom_barcode_counts = collections.Counter()

    for align in bam.fetch(contig=chrom):

        prefix_seq = align.get_forward_sequence()[: args.window]
        prefix_qv = align.get_forward_qualities()[: args.window]

        p_alignment = parasail_alg(
            s1=prefix_seq,
            s2=probe_seq,
            open=args.gap_open,
            extend=args.gap_extend,
            matrix=matrix,
        )

        adapter1_ed, barcode, bc_start = parse_probe_alignment(
            p_alignment, adapter1_probe_seq, args
        )

        # Require minimum read1 edit distance
        if adapter1_ed <= args.max_adapter1_ed:
            qscores = find_feature_qscores(barcode, p_alignment, prefix_seq, prefix_qv)
            chrom_barcode_counts[barcode] += 1

            # Strip out insertions from alignment to get read barcode sequence
            barcode = barcode.replace("-", "")

            # Uncorrected cell barcode = CR:Z
            align.set_tag("CR", barcode, value_type="Z")
            # Cell barcode quality score = CY:Z
            align.set_tag("CY", qscores, value_type="Z")

            # print(read_id)
            # print(p_alignment.traceback.ref)
            # print(p_alignment.traceback.comp)
            # print(p_alignment.traceback.query)
            # print()

            # Only write BAM entry in output file if it will have CR and CY tags
            bam_out.write(align)

    bam.close()
    bam_out.close()

    return chrom_bam.name, chrom_barcode_counts


def load_superlist(superlist):
    """
    Read contents of the file containing all possible cell barcode sequences.
    File can be uncompressed or gzipped.

    :param superlist: Path to file containing all possible cell barcodes, e.g.
        3M-february-2018.txt
    :type superlist: str
    :return: Set of all possible cell barcodes
    :rtype: set
    """
    ext = pathlib.Path(superlist).suffix
    fn = pathlib.Path(superlist).name
    wl = []
    if ext == ".gz":
        with gzip.open(superlist, "rt") as file:
            for line in tqdm(file, desc=f"Loading barcodes in {fn}", unit=" barcodes"):
                wl.append(line.strip())
    elif ext == ".txt":
        with open(superlist) as file:
            for line in tqdm(file, desc=f"Loading barcodes in {fn}", unit=" barcodes"):
                wl.append(line.strip())
    wl = set(wl)
    return wl


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

    # logger.info("Loading barcode superlist")
    wl = load_superlist(args.superlist)

    # Create temporary directory
    if os.path.exists(args.tempdir):
        shutil.rmtree(args.tempdir, ignore_errors=True)
    os.mkdir(args.tempdir)

    # Process BAM alignments from each chrom separately
    logger.info(f"Extracting uncorrected barcodes from {args.bam}")
    func_args = []
    for chrom in chroms.keys():
        func_args.append((args.bam, chrom, args))

    results = launch_pool(align_adapter, func_args, args.threads)
    chrom_bam_fns, chrom_barcode_counts = list(zip(*results))
    barcode_counts = sum(chrom_barcode_counts, collections.Counter())

    # Filter barcode counts against barcode superlist
    logger.info(f"Writing superlist-filtered barcode counts to {args.output_barcodes}")
    f_barcode_counts = open(args.output_barcodes, "w")
    barcode_counts_sorted = sorted(
        barcode_counts.items(), key=lambda item: item[1], reverse=True
    )
    for barcode, n in barcode_counts_sorted:
        if barcode in wl:
            f_barcode_counts.write(f"{barcode}\t{n}\n")
    f_barcode_counts.close()

    logger.info(f"Writing BAM with uncorrected barcode tags to {args.output_bam}")
    tmp_bam = tempfile.NamedTemporaryFile(
        prefix="tmp.align.", suffix=".unsorted.bam", dir=args.tempdir, delete=False
    )
    merge_parameters = ["-f", tmp_bam.name] + list(chrom_bam_fns)
    pysam.merge(*merge_parameters)

    pysam.sort("-@", str(args.threads), "-o", args.output_bam, tmp_bam.name)

    logger.info("Cleaning up temporary files")
    shutil.rmtree(args.tempdir, ignore_errors=True)


if __name__ == "__main__":
    args = parse_args()

    main(args)
