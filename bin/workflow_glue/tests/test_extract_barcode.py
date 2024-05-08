"""Test adapter_scan_vsearch."""
import functools
import os
import tempfile
from unittest.mock import Mock

import pandas as pd
import pytest
from pytest import fixture
from workflow_glue import extract_barcode
from ..sc_util import rev_cmp  # noqa: ABS101

# prevent stdout writing in align_adapter even with `pytest -s`
devnull = open(os.devnull, 'w')
extract_barcode.align_adapter = functools.partial(
    extract_barcode.align_adapter, fastq_out=devnull)


def gene():
    """Get a randomly-generated gene seuence."""
    return (
        "ATTCAGCGCTGGAGACCGAGCGCCCCGCAAAGGGCCTGATCT"
        "ATCGCGCACGGGACTACTCATTGGGACTGCGGCAATAGGGGAGGGGCCTAACAACGTT")


def make_fastq(
        read_adapter1='CTACACGACGCTCTTCCGATCT',
        read_barcode='AAACCCAAGAAACACT',
        read_umi='GACTGACTGACT',
        read_polyt='T'*12,
        rev=True):
    """Create a synthetic fastq file containing adapter and barcode sequences."""
    read = \
        f'{read_adapter1}{read_barcode}{read_umi}{read_polyt}{gene()}'
    if rev:
        read = rev_cmp(read)

    # Make a FASTQ file containing the read and a quality score of 60.
    fastq = (
        '@test_id\n'
        f"{read}\n"
        '+\n'
        f"{'?' * len(read)}\n"
    )

    # Write out a test FASTQ
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.fastq', delete=False) as fh_fq:
        fh_fq.write(fastq)
        fq_fn = fh_fq.name

    return fq_fn


@fixture
def make_superlist():
    """Make a small superlist(whitelist) of barcodes."""
    superlist = (
        "AAACCCAAGAAACACT\n"
        "AAACCCAAGAAACCAT")

    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv', delete=False) as fh_sl:
        fh_sl.write(superlist)
        superlist_fname = fh_sl.name
    return superlist_fname


@fixture()
def args(make_superlist):
    """Mock Args with workflow defaults set."""
    class Args:
        contig = 'chr17'
        fastq = make_fastq()
        match = 5
        mismatch = -1
        acg_to_n_match = 1
        t_to_n_match = 1
        adapter1_seq = 'CTACACGACGCTCTTCCGATCT'
        adapter1_suff_length = 10
        kit = '3prime'
        barcode_length = 16
        umi_length = 12
        window = 100
        gap_open = 2
        gap_extend = 4
        max_adapter1_ed = 3
        min_barcode_qv = 15
        polyt_length = 10
        superlist = make_superlist
        verbosity = 2

    return Args


def test_main(args):
    """Test the final output from main()."""
    # Make temp files to store the output
    counts_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')
    tags_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')
    trimmed_fastq_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')

    args.fastq = make_fastq(read_adapter1=args.adapter1_seq, rev=True)
    args.output_barcode_counts = counts_file.name
    args.output_read_tags = tags_file.name
    args.output_trimmed_fastq = trimmed_fastq_file.name
    extract_barcode.main(args)

    # Barcode we expect to find in the input BAM
    # For 3prime this will be reverse
    expected_barcode = 'AAACCCAAGAAACACT'

    counts_result = pd.read_csv(counts_file.name, sep='\t')
    assert counts_result.shape == (1, 2)
    assert counts_result.iat[0, 0] == expected_barcode
    assert counts_result.iat[0, 1] == 1

    tags_result = pd.read_csv(tags_file.name, sep='\t', index_col=0)

    # Check we have correct number of rows returned # read_id(index), CR, CY, UR, UY
    assert tags_result.shape == (1, 4)
    assert tags_result.loc['test_id', 'CR'] == 'AAACCCAAGAAACACT'
    assert tags_result.loc['test_id', 'UR'] == 'GACTGACTGACT'

    # TODO: test if barcode missing from superlist


@pytest.mark.parametrize(
    'adapter1_seq,tags_results_shape,counts_results_shape',
    [
        # Tags files should have 5 columns (read_id, CR, CY, UR, UY) ,
        # counts results should have one column (count) with read_id index
        ['CTACACGACGCTCTTCCGATCT', (1, 5), (1, 1)],  # ED 0
        ['CTACACGACGCTCTTCCGAggg', (1, 5), (1, 1)],  # ED 3
        ['CTACACGACGCTCTTCCGgggg', (0, 5), (0, 1)]   # ED 4; no results
    ]
)
def test_align_adapter(args, adapter1_seq, tags_results_shape, counts_results_shape):
    """Test the identification of adapter1 sequences.

    algin_adapter() should return results with a max adapter1 edit distance,
    which defaults to 3.
    """
    tags_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')
    trimmed_fastq_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')
    args.fastq = make_fastq(read_adapter1=adapter1_seq, rev=True)
    args.output_read_tags = tags_file.name
    args.kit = '3prime'
    args.output_trimmed_fastq = trimmed_fastq_file.name
    df_counts = extract_barcode.align_adapter(args)
    assert df_counts.shape == counts_results_shape

    df_tags = pd.read_csv(tags_file.name, sep='\t')
    assert df_tags.shape == tags_results_shape


def ascii_decode_qscores(string):
    """Convert ASCII character quality values into integers."""
    return list(map(lambda x: ord(x) - 33, string))


def ascii_encode_qscores(integers):
    """Convert integer quality values into ASCII characters."""
    return "".join(map(lambda x: chr(x + 33), integers))


@pytest.mark.parametrize(
    # 3 prime tests
    'query,query_aln,query_ascii_q,expected_adapter1_ed',
    [
        # Adapter 1             BC                UMI          polyT

        # 100% match of the query adapter1 and the 10bp adapter1 prefix in the ref probe
        ["CTACACGACGCTCTTCCGATCT AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "            CTTCCGATCT AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "?????????????????????? ()*+,-./01234567 89:;<=>?@ABC ????????????",
         0],

        # 2 bp substitution in the adapter1 query
        ["CTACACGACGCTCTTCCGATaa AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "            CTTCCGATaa AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "?????????????????????? ()*+,-./01234567 89:;<=>?@ABC ????????????",
         2],

        # 2 bp deletion in the adapter1 query
        ["CTACACGACGCTCTTCCGAT   AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "            CTTCCGAT-- AAACCCAAGAAACACT GACTGACTGACT TTTTTTTTTTTT",
         "????????????????????   ()*+,-./01234567 89:;<=>?@ABC ????????????",
         2],

    ]
)
def test_parse_probe_alignment(query, query_aln, query_ascii_q, expected_adapter1_ed):
    """Test_parse_probe_alignment.

    In this test, a mocked parasail alignment is created. We want to test that the
    correct barcode, UMI and associated quality scores are extracted from the
    query.

    :param query: read query
    :param query_aln: the query alignment that would result from parasail alignment to
        the reference probe
    :param query_ascii_q: the ascii-encoded qualitey string associated with the query
    :param: expected_adapter1_ed: the expected edit distance of the adapter1
    """
    # Build a mock parasail alignment result. Although there would be geen sequrnce
    # after the polyT, we can omit it here.

    # This is the read including the full 22bp adapter1 probe
    #       adapter1                BC               UMI          PolyT
    barcode, umi = query.split()[1:3]
    barcode_q, umi_q = query_ascii_q.split()[1:3]
    query = query.replace(' ', '')
    query_aln = query_aln.replace(' ', '')
    query_ascii_q = query_ascii_q.replace(' ', '')
    qual_ints = ascii_decode_qscores(query_ascii_q)

    # The parasail reference alignment. Contains only the 10 bp suffix of the adapter1
    # For 3prime kit
    ref_align = (
        # 10 bp A1  Ns for BC        Ns for UMI   PolyT
        "CTTCCGATCT NNNNNNNNNNNNNNNN NNNNNNNNNNNN TTTTTTTTTTTT"
    ).replace(' ', '')

    p_alignment = Mock()
    p_alignment.traceback.query = query_aln
    p_alignment.traceback.ref = ref_align

    adapter1_probe_suffix = 'CTTCCGATCT'  # Forward seq

    (
        adapter1_editdist, barcode_result, umi_result,
        bc_qscores, umi_qscores
    ) = extract_barcode.parse_probe_alignment(
        p_alignment, adapter1_probe_suffix,
        16, 12, qual_ints, query)

    # convert the return qscores to ascii-encoded values
    bc_qscores = ascii_encode_qscores(bc_qscores)
    umi_qscores = ascii_encode_qscores(umi_qscores)

    assert adapter1_editdist == expected_adapter1_ed
    assert barcode_result == barcode
    assert umi_result == umi
    assert bc_qscores == barcode_q
    assert umi_qscores == umi_q
