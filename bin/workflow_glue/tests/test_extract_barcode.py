"""Test adapter_scan_vsearch."""

import subprocess as sub
import tempfile
from unittest.mock import Mock

import pandas as pd
import pytest
from pytest import fixture
from workflow_glue.extract_barcode import (
    align_adapter, main, parse_probe_alignment)


def gene():
    """Get a randomly-generated gene seuence."""
    return (
        "ATTCAGCGCTGGAGACCGAGCGCCCCGCAAAGGGCCTGATCT"
        "ATCGCGCACGGGACTACTCATTGGGACTGCGGCAATAGGGGAGGGGCCTAACAACGTT")


def make_bam(
        read_adapter1='CTACACGACGCTCTTCCGATCT',
        read_barcode='AAACCCAAGAAACACT',
        read_umi='GACTGACTGACT',
        read_polyt='T'*12):
    """Create a synthetic bam file containing adapter and barcode sequences."""
    # make a bam
    header = """@SQ	SN:chr17	LN:10000000\n"""

    read = \
        f'{read_adapter1}{read_barcode}{read_umi}{read_polyt}{gene()}'

    # Make a sam file containing the read and a quality qscore of 60.
    sam = (
        f'{header}'
        f"test_id\t0\tchr17\t1\t60\t{len(read)}M\t*\t0\t0\t"
        f"{read}\t{'?' * len(read)}"
    )

    # Write out a test BAM
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.sam', delete=False) as fh_sam:
        fh_sam.write(sam)
        sam_file = fh_sam.name

    bam = 'test.bam'
    sub.check_output(['samtools', 'view', sam_file, '-o', bam])
    sub.check_output(['samtools', 'index', bam])

    return bam


@fixture
def make_superlist():
    """Make a samll superlist(whitelist) of barcodes."""
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
        bam = make_bam
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
    counts_file = tempfile.NamedTemporaryFile(
            mode='w', suffix='.tsv')
    tags_file = tempfile.NamedTemporaryFile(
        mode='w', suffix='.tsv')

    args.bam = make_bam(read_adapter1=args.adapter1_seq)

    args.output_barcode_counts = counts_file.name
    args.output_read_tags = tags_file.name
    main(args)

    expected_barcode = 'AAACCCAAGAAACACT'

    counts_result = pd.read_csv(
        counts_file.name, sep='\t', names=['barcode', 'count'])
    assert counts_result.shape == (1, 2)
    assert counts_result.iat[0, 0] == expected_barcode
    assert counts_result.iat[0, 1] == 1

    tags_result = pd.read_csv(tags_file.name, sep='\t', index_col=0)

    assert tags_result.shape == (1, 8)
    assert tags_result.loc['test_id', 'CR'] == 'AAACCCAAGAAACACT'
    assert tags_result.loc['test_id', 'UR'] == 'GACTGACTGACT'

    # TODO: test if barcode missing from superlist


@pytest.mark.parametrize(
    'adapter1_seq,tags_results_shape,counts_results_shape',
    [
        ['CTACACGACGCTCTTCCGATCT', (1, 9), (1, 1)],  # ED 0
        ['CTACACGACGCTCTTCCGAggg', (1, 9), (1, 1)],  # ED 3
        ['CTACACGACGCTCTTCCGgggg', (0, 9), (0, 1)]   # ED 4
    ]
)
def test_align_adapter(args, adapter1_seq, tags_results_shape, counts_results_shape):
    """Test the identification of adapter1 sequences.

    algin_adapter() should return results with a max adapter1 edit distance,
    which defaults to 3.

    How do we change adapter1 seq?
    """
    args.bam = make_bam(read_adapter1=adapter1_seq)
    df_tags, df_counts = align_adapter(args)
    assert df_tags.shape == tags_results_shape
    assert df_counts.shape == counts_results_shape


def ascii_encode_qscores(integers):
    """Convert integer quality values into ASCII characters."""
    return "".join(map(lambda x: chr(x + 33), integers))


@pytest.mark.parametrize(
    'align_adapter1,align_barcode,align_umi,adapter1_edit_dist',
    [
        # Full adpt1 match, no deletions in BC/UMI
        ['CTTCCGATCT', 'AAACCCAAGAAACACT', 'GACTGACTGACT', 0],
        # 3 ED for adpt1
        ['gaaCCGATCT', 'AAACCCAAGAAACACT', 'GACTGACTGACT', 3],
        # Deletion in barcode
        ['CTTCCGATCT', 'AAACCCAAGAAACA-T', 'GACTGACTGACT', 0],
        # Deletion in UMI
        ['CTTCCGATCT', 'AAACCCAAGAAACA-T', '-ACTGACTGACT', 0],
        # Deletion in barcode and umi
        ['CTTCCGATCT', 'A--CCCAAGAAACACT', 'GA--GACTGACT', 0],
        # Deletion in adapter, barcode, and umi
        ['CTTCCGAT-T', 'A--CCCAAGAAACACT', 'GA--GACTGACT', 1],

    ]

)
def test_parse_probe_alignment(
        args, align_adapter1, align_barcode, align_umi, adapter1_edit_dist):
    """Test_parse_probe_alignment.

    Mock various parasail alignments to test alignment parsing.
    """
    # Some defaults for the 3prime kit. These are passed to parse_probe_alignemtn
    barcode_length = 16
    umi_length = 12

    # Get real barcode and UMI by stripping any deletion characters form the alignment
    actual_bc = align_barcode.replace('-', '')
    actual_bc_len = len(actual_bc)
    actual_umi = align_umi.replace('-', '')
    actual_umi_len = len(actual_umi)
    actual_a1_len = len(align_adapter1.replace('-', ''))

    # build probe - only the 10bp suffix of the adapter1 is included in the probe.
    adapter1_probe_suffix = 'CTTCCGATCT'

    p_alignment = Mock()
    reference_probe = \
        f"{adapter1_probe_suffix}{'N'* barcode_length}{'N'* umi_length}{'T' * 12}"
    p_alignment.traceback.ref = reference_probe

    # Build query alignment from the first 100bp of the read
    read_prefix = f"{align_adapter1}{align_barcode}{align_umi}{'T' * 20}{gene}"[:100]
    # The aligned query sequence
    p_alignment.traceback.query = (
        f"{align_adapter1}{align_barcode}{align_umi}{'T' * 12}")
    # The read prefix used as query
    p_alignment.query = read_prefix
    # The last aligned position in the query
    p_alignment.end_query = \
        len(align_adapter1) + actual_bc_len + actual_umi_len + 12 - 1

    aln_quality = range(len(read_prefix))

    actual_bc_qual = ascii_encode_qscores(
        aln_quality[actual_a1_len: actual_a1_len + actual_bc_len])
    actual_umi_qual = ascii_encode_qscores(aln_quality[
        actual_a1_len + actual_bc_len: actual_a1_len + actual_bc_len + actual_umi_len])

    (
        adapter1_editdist, barcode_result, umi_result,
        bc_qscores, umi_qscores
    ) = parse_probe_alignment(
            p_alignment, adapter1_probe_suffix,
            barcode_length, umi_length, aln_quality)

    assert adapter1_editdist == adapter1_edit_dist
    assert barcode_result == actual_bc
    assert umi_result == actual_umi
    assert bc_qscores == actual_bc_qual
    assert umi_qscores == actual_umi_qual
