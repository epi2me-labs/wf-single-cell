"""Adapter scan vsearch."""
import collections
from pathlib import Path
import subprocess
import sys

import pandas as pd
import polars as pl
import pysam

from .sc_util import kit_adapters, StatsSummary  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101

logger = get_named_logger("AdaptScan")

# adapter information
compat_adapters = {"adapter1_f": "adapter2_f", "adapter2_r": "adapter1_r"}
configs = {
    "adapter2_r-adapter2_f": 'double_adapter2',
    "adapter2_f-adapter2_r": 'double_adapter2',
    "adapter1_r-adapter1_f": 'double_adapter1',
    "adapter1_f-adapter1_r": 'double_adapter1',
    "adapter2_f": 'single_adapter2',
    "adapter2_r": 'single_adapter2',
    "adapter1_f": 'single_adapter1',
    "adapter1_r": 'single_adapter1',
    "*": "no_adapters"}

# for reverse complementation
complement_trans = str.maketrans("ACGTacgt", "TGCAtgca")


def argparser():
    """Create argument parser."""
    parser = wf_parser("adapt_scan")
    parser.add_argument(
        "fastq", type=Path,
        help="FASTQ of ONT reads.")

    parser.add_argument(
        "--per_read_summary", type=Path,
        help="Output file name for per-read adapter configurations.")

    parser.add_argument(
        "--summary", type=Path,
        help="Output for JSON summary of adapter classifications.")

    parser.add_argument(
        "--kit",
        help="Specify either the 10X 3' gene expression kit (3prime), the 5' \
        gene expression kit (5prime), or the multiome kit (multiome) This \
        determines which adapter sequences to search for in the reads.",
        default="3prime", choices=['3prime', '5prime', 'multiome', 'visium'])

    parser.add_argument(
        "--min_adapter_id", type=float, default=0.7,
        help="Minimum adapter alignment identity for VSEARCH.")

    parser.add_argument(
        "--keep_fl_only",
        help="Only write full length reads.",
        action="store_true", default=False)

    parser.add_argument(
        "--adapters_fasta", type=Path, default="adapter_seqs.fasta",
        help="Filename for adapter query sequences.")

    parser.add_argument(
        "--threads", type=int, default=8,
        help="Compute threads for vsearch.")
    return parser


def write_adapters_fasta(adapter1_seq, adapter2_seq, output):
    """Write adapters fasta for use with vsearch."""
    logger.info(f"Writing adapter sequences to {output}.")
    adapters = {
        "adapter1_f": adapter1_seq,
        "adapter1_r": adapter1_seq[::-1].translate(complement_trans),
        "adapter2_f": adapter2_seq,
        "adapter2_r": adapter2_seq[::-1].translate(complement_trans)}
    with open(output, 'w') as fh:
        for name, seq in adapters.items():
            fh.write(f">{name}\n{seq}\n")


def call_vsearch(fastq, output, min_adapter_id, adapters_fasta, threads):
    """Call vsearch."""
    vsearch_cmd = f"seqkit fq2fa {fastq} | vsearch --usearch_global - \
        --db {adapters_fasta} --minseqlength 20 --maxaccepts 5 --id {min_adapter_id} \
        --strand plus --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
        'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
        --userout {output} --threads {threads}"

    logger.info(vsearch_cmd)
    p = subprocess.Popen(
        vsearch_cmd, shell=True, stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    p.communicate()
    return None


def get_valid_adapter_pair_positions_in_read(hits):
    """Get valid adapter positions.

    Note on vsearch terminology
    qilo:  the first query position aligning with the reference
    qihi: the last query position aligning with the reference
    vsearch results.

    :param hits: vsearch results
    :type hits: polars.DataFrame
    """
    read_id = hits.item(0, 'query')

    def _create_pair(
            first_adapter, second_adapter, adapter_1_idx, adapter_2_idx, pair_id):
        # Get region bounded by adapter pair
        start = hits.item(adapter_1_idx, 'qilo') - 1
        end = hits.item(adapter_2_idx, 'qihi')

        # Adapter2 is not needed for downstream processing
        # So trim adapter2_F or adapter2_R but keep adapter_1
        if first_adapter == 'adapter1_f':
            # Adapter2 is at the end of the read, set end to be just before it
            end = hits.item(adapter_2_idx, 'qilo') - 1
        elif first_adapter == 'adapter2_r':
            start = hits.item(adapter_1_idx, 'qihi')
        else:
            raise ValueError("Unexpected adapter pairing produced")

        fl_pair = {
            "read_id": f"{read_id}_{pair_id}",
            "config": f"{first_adapter}-{second_adapter}",
            "start": start,
            "end": end}
        fl_pair["strand"] = "+" if first_adapter == "adapter1_f" else "-"
        return fl_pair

    # the polars row searches are slow in simple cases
    if len(hits) == 1:
        return list()
    elif len(hits) == 2:
        for first, second in compat_adapters.items():
            if hits["target"][0] == first and hits["target"][1] == second:
                pair = _create_pair(first, second, 0, 1, 0)
                return [pair]
        return list()  # there was no valid pair
    else:
        # traverse the hits looking for pairs
        pair_id = 0
        fl_pairs = list()
        for first, second in compat_adapters.items():
            adapter_1_idxs = hits.select(
                [pl.arg_where(pl.col("target") == first)]).to_series()
            for adapter_1_idx in adapter_1_idxs:
                # check the next hit, is it the corresponding adapter from the pair
                adapter_2_idx = adapter_1_idx + 1
                if adapter_2_idx < len(hits):
                    if hits.item(adapter_2_idx, 'target') == second:
                        fl_pairs.append(
                            _create_pair(
                                first, second, adapter_1_idx, adapter_2_idx, pair_id))
                        pair_id += 1
        return fl_pairs
    raise RuntimeError("Should not reach here")


def parse_vsearch(tmp_vsearch):
    """Parse vsearch adapter scan results.

    Parse the vsearch results, identifying adapter configurations including
    pairs of compatible adapters which represent full length read segments.

    :param tmp_vsearch: TSV-formatted vsearch results from 10x adapter search.
    :type tmp_vsearch: str
    :return: read data including adapter configurations, strand and
        full length read segment locations.
    :rtype: dict
    """
    logger.info("Reading data")
    vsearch_colnames = [
        "query", "target", "id", "alnlen", "mism", "opens",
        "qilo", "qihi", "qstrand", "tilo", "tihi",
        "ql", "tl"]

    schema = {
        'query': pl.Utf8,
        'target': pl.Categorical,
        'qilo': pl.UInt32,
        'qihi': pl.UInt32,
        'ql': pl.UInt32}

    columns_idxs = [vsearch_colnames.index(x) for x in schema.keys()]
    df = pl.read_csv(
        source=tmp_vsearch,
        has_header=False,
        separator='\t',
        dtypes=schema,
        columns=columns_idxs,
        new_columns=list(schema.keys())
    ).sort(["query", "qilo"])
    logger.info("Finished reading and sorting data")

    # grouby isn't terribly efficient even when sorted,
    # we can do better with simple linear chunking
    def _yield_reads():
        last_i, last = 0, df["query"][0]
        for i, name in enumerate(df["query"][1:], start=1):
            if name != last:
                yield last, df[last_i:i]
                last_i, last = i, name
        yield last, df[last_i:]

    # find adapters within reads by grouping hits to single read
    read_info = collections.defaultdict(list)  # key is original read ID
    for orig_read_id, read_result in _yield_reads():
        orig_adapter_config = "-".join(read_result["target"])

        # find valid consecutive adapter pairs, this allows for some
        # forms of concatemeric reads to be recovered
        fl_pairs = get_valid_adapter_pair_positions_in_read(read_result)

        if len(fl_pairs) > 0:
            # valid pairs found
            for fl_pair in fl_pairs:
                read_info[orig_read_id].append({
                    "orig_read_id": orig_read_id,
                    "read_id": fl_pair["read_id"],
                    "readlen": fl_pair["end"] - fl_pair["start"],
                    "start": fl_pair["start"],
                    "end": fl_pair["end"],
                    "fl": True,
                    "stranded": True,
                    "orig_strand": fl_pair["strand"],
                    "orig_adapter_config": orig_adapter_config,
                    "adapter_config": fl_pair["config"],
                    "lab": "full_len"})
        else:
            # no valid adapter pairs found, try to recover reads that have
            # one useful adapter it.
            read_row = read_result.row(0, named=True)
            read_id = f"{orig_read_id}_0"
            fl = False
            stranded = False
            strand = "*"
            readlen = read_row['ql']
            start = 0
            end = readlen
            adapter_config = "-".join(read_result["target"])
            config_type = configs.get(adapter_config, 'other')

            # If there is only a single adapter, we can assign a strand and
            # also trim the adapter side of the read. If there are any other
            # configurations, such as multiple unpaired adapters, we do not
            # attempt to process further.
            if config_type == 'single_adapter2':
                # We want to strand and trim reads where we only have an
                # adapter2 sequence but no adapter1. These MIGHT contain the
                # cell barcode, UMI, polyT (for --kit 3prime) and cDNA
                # sequence, but since the adapter1 is low-quality, these might
                # be of dubious value. We can only trim the adapter2 end of the
                # read based on adapter2 aligned positions, but won't touch the
                # putative adapter1 end.
                stranded = True
                if adapter_config == "adapter2_f":
                    strand = "+"
                    start = 0
                    end = read_row['qihi']
                    readlen = end - start
                elif adapter_config == "adapter2_r":
                    strand = "-"
                    start = read_row['qilo'] - 1
                    end = read_row['ql']
                    readlen = end - start
                else:
                    raise Exception("Unexpected adapter config!")
            elif config_type == "single_adapter1":
                # We want to strand and trim reads where we only have a
                # adapter1 sequence but no adapter2. These should contain the
                # cell barcode, UMI, polyT (if --kit 3prime) and cDNA sequence,
                # so should still have significant value. We can only trim the
                # adapter1 end of the read but won't touch the putative
                # adapter2 end.
                stranded = True
                if adapter_config == "adapter1_f":
                    strand = "+"
                    start = read_row['qilo'] - 1
                    end = read_row['ql']
                    readlen = end - start
                elif adapter_config == "adapter1_r":
                    strand = "-"
                    start = 0
                    end = read_row['qihi']
                    readlen = end - start
                else:
                    raise Exception("Unexpected adapter config!")

            read_info[orig_read_id].append({
                "orig_read_id": orig_read_id,
                "read_id": read_id,
                "readlen": readlen,
                "start": start,
                "end": end,
                "fl": fl,
                "stranded": stranded,
                "orig_strand": strand,
                "orig_adapter_config": orig_adapter_config,
                "adapter_config": adapter_config,
                "lab": config_type})

    return read_info


def make_read_table(read_info):
    """Make per-read data table."""
    return pd.DataFrame(read_info).set_index('readlen', drop=True)


class AdapterSummary(StatsSummary):
    """Summary dictionary for storing adapter configuration summaries."""

    fields = {
        "subreads", "full_length", "stranded", "plus", "minus",
        "single_adapter1", "double_adapter1",
        "single_adapter2", "double_adapter2",
        "no_adapters", "other"}

    @classmethod
    def from_pandas(cls, df):
        """Create an instance from a pandas dataframe."""
        stats = dict()
        stats["subreads"] = len(df)
        stats["full_length"] = len(df[df.fl].index)
        stats["stranded"] = len(df[df['orig_strand'] != '*'])
        stats["plus"] = len(df[df['orig_strand'] == '+'])
        stats["minus"] = len(df[df['orig_strand'] == '-'])
        # summary config: "single_adapter1" etc
        stats.update(df['lab'].value_counts().astype(int).to_dict())
        # don't worry about precise configuration
        return cls(stats)


def create_stranded_reads(fastq, read_info, kit, fl_only):
    """Write stranded and trimmed fastq.

    For reads with adapter1, trim up to the adapter. The adapter sequence is needed
    in the read for downstream BC searching.

    For reads containing an adapter2 sequence, trim that away as it's not needed.
    """
    adap_rev_config = {
        "adapter1_f": "adapter1_r",
        "adapter1_r": "adapter1_f",
        "adapter2_f": "adapter2_r",
        "adapter2_r": "adapter2_f"}

    # Iterate through FASTQ reads and re-write them with proper stranding based on
    # the results of the VSEARCH alignments.
    with pysam.FastxFile(fastq, persist=False) as f_in:
        for entry in f_in:
            if read_info.get(entry.name):
                # This read had some VSEARCH hits for adapter sequences
                for subread in read_info[entry.name]:
                    if fl_only and not subread["fl"]:
                        continue
                    # The read may contain more than one subread -
                    # segments flaked by compatible adapters. Extract these and
                    # trim
                    subread_id = subread["read_id"]
                    subread_seq = entry.sequence[subread["start"]: subread["end"]]
                    subread_quals = entry.quality[subread["start"]: subread["end"]]

                    if any([
                        (subread["orig_strand"] == "-" and kit == '5prime'),
                        (subread["orig_strand"] == '+' and kit in [
                            '3prime', 'multiome', 'visium'])
                    ]):
                        # Do any necessary stranding
                        #
                        # 5prime read structure:
                        #   adapter1-BC-UMI-TSO-cDNA-polyA
                        # 5 prime kit adapters are the same sense as the
                        # mRNA, therefore the reads are re-oriented if the original
                        # strand is '-'
                        #
                        # 3prime read structure:
                        #   adapter1-BC-UMI-polyT-cDNA
                        # For 3prime and multiome kits, the adapters are on the
                        # opposite strand relative to mRNA sense, so reverse
                        # complement the read if the original strand is '+'
                        #
                        # Reverse the config string and reverse
                        # complement the subread
                        rc_config = "-".join(
                            [adap_rev_config[adapter] for adapter in
                             subread["adapter_config"].split("-")[::-1]])

                        subread["adapter_config"] = rc_config
                        subread_seq = subread_seq[::-1].translate(
                            complement_trans)
                        subread_quals = subread_quals[::-1]

                    if not subread_seq:
                        # Rarely a sequence is totally trimmed away.
                        continue
                    yield f"@{subread_id}\n{subread_seq}\n+\n{subread_quals}\n"


def main(args):
    """Entry point."""
    # TODO: this is a little weird \:D/
    adapters = kit_adapters[args.kit]
    args.adapter1_seq = adapters['adapter1']
    args.adapter2_seq = adapters['adapter2']
    adapter_file = 'adapters.fasta'
    write_adapters_fasta(
        adapters['adapter1'], adapters['adapter2'], adapter_file)

    logger.info("Running vsearch")
    vsearch_out = args.fastq.with_suffix(".vsearch.tsv")
    call_vsearch(
        args.fastq, vsearch_out, args.min_adapter_id, adapter_file, args.threads)

    logger.info("Parsing vsearch hits.")
    read_info = parse_vsearch(vsearch_out)

    logger.info(f"Creating fastq for {len(read_info)} reads.")
    read_data = create_stranded_reads(
        args.fastq, read_info, args.kit, args.keep_fl_only)
    for read in read_data:
        sys.stdout.write(read)

    if args.per_read_summary is not None or args.summary is not None:
        logger.info("Making summary documents.")
        data = list()
        for subreads in read_info.values():
            data.extend(subreads)
        df = pd.DataFrame(data).set_index('read_id', drop=True)
        if args.per_read_summary:
            logger.info(f"Writing output table to {args.per_read_summary}")
            df.to_csv(args.per_read_summary, sep="\t", index=True)
        # TODO: make the summary without going via table, this could be done directly
        #       in parse_vsearch()
        if args.summary:
            logger.info(f"Writing JSON summary to {args.summary}")
            summary = AdapterSummary.from_pandas(df)
            summary.to_json(args.summary)

    logger.info("Finished.")
