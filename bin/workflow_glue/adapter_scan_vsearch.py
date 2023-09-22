#!/usr/bin/env python
"""Adapter scan vsearch."""
import gzip
from pathlib import Path
import subprocess
import tempfile

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import polars as pl
import pysam
from workflow_glue.sc_util import kit_adapters

from .util import get_named_logger, wf_parser  # noqa: ABS101


compat_adapters = {"adapter1_f": "adapter2_f", "adapter2_r": "adapter1_r"}

# Map adapter configurations to configuration class.
configs = {
    "adapter2_r-adapter2_f": 'double_adapter2',
    "adapter2_f-adapter2_r": 'double_adapter2',
    "adapter1_r-adapter1_f": 'double_adapter1',
    "adapter1_f-adapter1_r": 'double_adapter1',
    "adapter2_f": 'single_adapter2',
    "adapter2_r": 'single_adapter2',
    "adapter1_f": 'single_adapter1',
    "adapter1_r": 'single_adapter1',
    "*": "no_adapters"
}


def argparser():
    """Create argument parser."""
    parser = wf_parser("adapt_scan")
    # Positional mandatory arguments
    parser.add_argument("fastq", help="FASTQ of ONT reads", type=Path)

    # Optional arguments
    parser.add_argument(
        "--output_fastq",
        help="Output file name for (gzipped) stranded FASTQ entries",
        type=Path,
        default="stranded.fastq.gz",
    )

    parser.add_argument(
        "--output_tsv",
        help="Output file name for adapter configurations",
        type=Path,
        default="adapters.tsv",
    )

    parser.add_argument(
        "-k",
        "--kit",
        help="Specify either the 10X 3' gene expression kit (3prime), the 5' \
        gene expression kit (5prime), or the multiome kit (multiome) This \
        determines which adapter sequences to search for in the reads",
        default="3prime",
        choices=['3prime', '5prime', 'multiome']
    )

    parser.add_argument(
        "-i",
        "--min_adapter_id",
        help="Minimum adapter alignment identity for VSEARCH",
        type=float,
        default=0.7,
    )

    parser.add_argument(
        "--only_strand_full_length",
        help="Do not try to strand-orient reads where either \
                        just a single adapter was found",
        action="store_true",
        default=False,
    )

    parser.add_argument(
        "-a",
        "--adapters_fasta",
        help="Filename for adapter query sequences",
        type=Path,
        default="adapter_seqs.fasta",
    )

    return parser


# complement translation table with support for regex punctuation
complement_trans = str.maketrans(
    "ACGTWSMKRYBDHVNacgtwsmkrybdhvn", "TGCAWSKMYRVHDBNtgcawskmyrvhdbn"
)

vsearch_colnames = [
    "query",
    "target",
    "id",
    "alnlen",
    "mism",
    "opens",
    "qilo",
    "qihi",
    "qstrand",
    "tilo",
    "tihi",
    "ql",
    "tl",
]


def write_adapters_fasta(adapter1_seq, adapter2_seq, output):
    """Write adapters fasta for use with vsearch."""
    adapters = []
    for adapter, seq in {
        "adapter1_f": adapter1_seq,
        "adapter1_r": adapter1_seq[::-1].translate(complement_trans),
        "adapter2_f": adapter2_seq,
        "adapter2_r": adapter2_seq[::-1].translate(complement_trans),
    }.items():
        entry = SeqRecord(Seq(seq), id=adapter, name="", description="")

        adapters.append(entry)
        with open(output, 'w') as fh:
            SeqIO.write(adapters, fh, "fasta")


def call_vsearch(fastq, min_adapter_id, adapters_fasta):
    """Call vsearch."""
    tmp_vsearch = fastq.with_suffix(".vsearch.tsv")

    vsearch_cmd = "seqkit fq2fa {fastq} | vsearch --usearch_global -  \
    --db {adapters} --minseqlength 20 --maxaccepts 5 --id {id} --strand plus \
    --wordlength 3 --minwordmatches 10 --output_no_hits --userfields \
    'query+target+id+alnlen+mism+opens+qilo+qihi+qstrand+tilo+tihi+ql+tl' \
    --userout {output}".format(
        fastq=fastq,
        id=min_adapter_id,
        adapters=adapters_fasta,
        output=tmp_vsearch,
    )

    p = subprocess.Popen(
        vsearch_cmd, shell=True, stderr=subprocess.PIPE,
        stdout=subprocess.PIPE)
    p.communicate()
    return tmp_vsearch


def get_valid_adapter_pair_positions_in_read(vs_result):
    """Get valid adapter positions.

    :param vs_result: vsearch results
    :type vs_result: polars.DataFrame
    """
    valid_pairs_n = 0
    fl_pairs = []

    # Find the first adapter of each segment.
    # If + strand, first adapter is adapter1_f
    # If - strand, first adapter is adapter2_r
    for first_adapter in compat_adapters.keys():
        adapter_1_idxs = vs_result.with_row_count().filter(
            vs_result["target"] == first_adapter)['row_nr']
        for adapter_1_idx in adapter_1_idxs:
            # For each found first adapter, examine next found adapter
            adapter_2_idx = adapter_1_idx + 1
            # Make sure there are enough alignments to allow this indexing
            if adapter_2_idx < vs_result.height:
                # Is the next found adapter an adapter2_f?
                if vs_result.item(adapter_2_idx, 'target') \
                        == compat_adapters[first_adapter]:
                    # This is a valid adapter pairing (adapter1_f-adapter2_f)
                    read_id = vs_result.item(0, 'query')
                    pair_str = (
                        f"{vs_result.item(adapter_1_idx, 'target')}-"
                        f"{vs_result.item(adapter_2_idx, 'target')}"
                    )
                    fl_pair = {
                        "read_id": f"{read_id}_{valid_pairs_n}",
                        "config": pair_str,
                        "start": vs_result.item(adapter_1_idx, 'qilo') - 1,
                        "end": vs_result.item(adapter_2_idx, 'qihi'),
                    }

                    fl_pair["strand"] = "+" if first_adapter == "adapter1_f" else "-"
                    valid_pairs_n += 1
                    fl_pairs.append(fl_pair)
    return fl_pairs


def parse_vsearch(tmp_vsearch, only_strand_full_length):
    """Parse vsearch adapter scan results.

    Parse the vsearch results, identifying adapter configurations including
    pairs of compatible adapters which represent full length read segments.

    :param tmp_vsearch: TSV-formatted vsearch results from 10x adapter search.
    :type tmp_vsearch: str
    :return: read data including adapter configurations, strand and
        full length read segment locations.
    :rtype: dict
    """
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
    )

    read_info = {}

    df = df.sort(["query", "qilo"])

    # Group the vsearch results by read_id. The resulting dataframe may
    # contain zero or more pairs of adapter hits representing full length
    # read segments. There may also be unpaired adapter hits where an
    # identified adapter does not have a consecutive compatible adapter.
    for orig_read_id, read_result in df.groupby("query"):
        # Get the string representation of the adapters identified in the read.
        # eg: 2 compatible forward adapters: 'adapter1_f-adapter2_f'
        orig_adapter_config = "-".join(read_result["target"])
        read_info[orig_read_id] = {}

        # Get any full length read segments by locating pairs of valid
        # consecutive adapters.
        # For example, an adapter1_f followed immediately by
        # an adapter2_f, or an adapter2_r followed immediately by
        # an adapter1_r.
        fl_pairs = get_valid_adapter_pair_positions_in_read(read_result)

        if len(fl_pairs) > 0:
            # There is at least one pair of compatible adapters
            # (a full length read segment). Add these to the read_info.
            for fl_pair in fl_pairs:
                read_info[orig_read_id][fl_pair['read_id']] = {
                    "readlen": fl_pair["end"] - fl_pair["start"],
                    "read_id": fl_pair["read_id"],
                    "start": fl_pair["start"],
                    "end": fl_pair["end"],
                    "fl": True,
                    "stranded": True,
                    "orig_strand": fl_pair["strand"],
                    "orig_adapter_config": orig_adapter_config,
                    "adapter_config": fl_pair["config"],
                    "lab": "full_len"}
        else:
            # No valid adapter pairs found. Either single adapter,
            # or weird artifact read.

            # The first vsearch result contains the info we need.
            read_row = read_result.row(0, named=True)
            # Make a single subread id.
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
                if not only_strand_full_length:
                    # We want to strand and trim reads where we only have
                    # an adapter2 sequence but no adapter1. These MIGHT contain
                    # the cell barcode, UMI, polyT (for --kit 3prime) and cDNA
                    # sequence, but since the adapter1 is low-quality, these
                    # might be of dubious value. We can only trim the
                    # adapter2 end of the vs_read_results based on adapter2
                    # aligned
                    # positions, but won't touch the putative adapter1 end.
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
                        raise Exception("Shouldn't be here!")
            elif config_type == "single_adapter1":
                if not only_strand_full_length:
                    # We want to strand and trim reads where we only have
                    # a adapter1 sequence but no adapter2. These should contain
                    # the cell barcode, UMI, polyT (if --kit 3prime) and cDNA
                    # sequence, so should still have significant value. We can
                    # only trim the adapter1 end of the vs_read_results,
                    # but won't touch the putative adapter2 end.
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
                        raise Exception("Shouldn't be here!")

            read_info[orig_read_id][read_id] = {
                "readlen": readlen,
                "read_id": read_id,
                "start": start,
                "end": end,
                "fl": fl,
                "stranded": stranded,
                "orig_strand": strand,
                "orig_adapter_config": orig_adapter_config,
                "adapter_config": adapter_config,
                "lab": config_type}

    return read_info


def write_tables(read_info, output_tsv):
    """Write the final table."""
    entries = []
    for _, entry in read_info.items():
        for _, sub_entry in entry.items():
            entries.append(sub_entry)
    df_table = pd.DataFrame(
        entries).set_index('readlen', drop=True)
    df_table.to_csv(output_tsv, sep="\t", index=True)


def write_stranded_fastq(fastq, read_info, output_fastq):
    """Write stranded fastq.

    Trimmed sub-read segments defined in read_info are written out with the
    correct orientation.
    """
    adap_rev_config = {
        "adapter1_f": "adapter1_r",
        "adapter1_r": "adapter1_f",
        "adapter2_f": "adapter2_r",
        "adapter2_r": "adapter2_f"
    }

    # Iterate through FASTQ reads and re-write them with proper stranding based
    # the results of the VSEARCH alignments.
    with pysam.FastxFile(fastq) as f_in:
        # with open(tmp_stranded_fastq, "w") as f_out:
        with gzip.open(output_fastq, "wb") as f_out:
            for entry in f_in:
                if read_info.get(entry.name):
                    # This read had some VSEARCH hits for adapter sequences
                    for subread_id, d in read_info[entry.name].items():
                        # The read may contain more than one subread -
                        # segments flaked by compatible adapters.
                        subread_seq = entry.sequence[d["start"]: d["end"]]
                        subread_quals = entry.quality[d["start"]: d["end"]]

                        if d["orig_strand"] == "-":
                            # Read in reverse orientation relatve to mRNA
                            # Reverse the config string and reverse
                            # complement the subread
                            rc_config = "-".join(
                                [adap_rev_config[a] for a in
                                 d["adapter_config"].split("-")[::-1]])

                            d["adapter_config"] = rc_config
                            subread_seq = subread_seq[::-1].translate(
                                complement_trans)
                            subread_quals = subread_quals[::-1]

                        f_out.write(
                            (f"@{subread_id}\n"
                             f"{subread_seq}\n"
                             f"+\n"
                             f"{subread_quals}\n").encode())


def main(args):
    """Entry point."""
    logging = get_named_logger("adapt_scan")
    adapters = kit_adapters[args.kit]
    args.adapter1_seq = adapters['adapter1']
    args.adapter2_seq = adapters['adapter2']

    # Create temp dir and add that to the args object
    p = Path(args.output_tsv)
    tempdir = tempfile.TemporaryDirectory(prefix="tmp.", dir=p.parents[0])
    args.tempdir = tempdir.name
    adapter_file = 'adapters.fasta'
    write_adapters_fasta(
        adapters['adapter1'], adapters['adapter2'], adapter_file)
    vsearch_results = call_vsearch(
        args.fastq, args.min_adapter_id, adapter_file)
    read_info = parse_vsearch(vsearch_results, args.only_strand_full_length)
    write_stranded_fastq(args.fastq, read_info, args.output_fastq)
    write_tables(read_info, args.output_tsv)
    logging.debug(f"Writing output table to {args.output_tsv}")


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
