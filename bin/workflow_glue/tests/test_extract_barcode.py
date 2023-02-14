"""Test adapter_scan_vsearch."""

import subprocess as sub
import tempfile
from unittest.mock import Mock

import numpy as np
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

    counts_result = pd.read_csv(
        counts_file.name, sep='\t', names=['barcode', 'count'])
    assert counts_result.shape == (1, 2)
    assert counts_result.iat[0, 0] == 'AAACCCAAGAAACACT'
    assert counts_result.iat[0, 1] == 1

    tags_result = pd.read_csv(tags_file.name, sep='\t', index_col=0)

    assert tags_result.shape == (1, 7)
    assert tags_result.loc['test_id', 'CR'] == 'AAACCCAAGAAACACT'
    assert tags_result.loc['test_id', 'UR'] == 'GACTGACTGACT'

    # TODO: test if barcode missing from superlist


@pytest.mark.parametrize(
    'adapter1_seq,tags_results_shape,counts_results_shape',
    [
        ['CTACACGACGCTCTTCCGATCT', (1, 8), (1, 1)],  # ED 0
        ['CTACACGACGCTCTTCCGAggg', (1, 8), (1, 1)],  # ED 3
        ['CTACACGACGCTCTTCCGgggg', (0, 8), (0, 1)]   # ED 4
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


@pytest.mark.parametrize(
    'align_adapter1,align_barcode,adapter1_edit_dist',
    [
        ['CTTCCGATCT', 'AAACCCAAGAAACACT', 0],  # full match
        ['CTTCCGAggg', 'AAACCCAAGAAACACT', 3],
        # ['CTTCCGATCT', 'AAACCCAAGAAA-ACT', 0],
    ]
)
def test_parse_probe_alignment(
        args, align_adapter1, align_barcode, adapter1_edit_dist):
    """Test_parse_probe_alignment.

    Mock various parasail alignments to test alignment parsing.
    """
    # Some defaults for the 3prime kit. These are passed to parse_probe_alignemtn
    barcode_length = args.barcode_length
    umi_length = args.umi_length
    umi = 'GACTGACTGACT'  # Use same for probe and alignment
    # Only the 10bp suffix of the adapter1 is included in the probe.
    adapter1_probe_suffix = 'CTTCCGATCT'

    # These go into making up the probe sequence and are condtant
    barcode_probe = 'AAACCCAAGAAACACT'
    umi_probe = 'GACTGACTGACT'
    polyt = 't' * 12  # same

    # Give each component of the read a unique quality score.
    # Only the 10bp adapter1 suffix is used in the probe
    read_adapter1_q = np.repeat(23, len(adapter1_probe_suffix))
    read_barcode_q = np.repeat(24, len(barcode_probe))
    read_umi_q = np.repeat(25, len(umi_probe))
    read_polyt_q = np.repeat(26, len(polyt))
    gene_q = np.repeat(27, len(gene()))

    prefix_qual = np.concatenate((
        read_adapter1_q,
        read_barcode_q,
        read_umi_q,
        read_polyt_q,
        gene_q
    ))

    read_umi_q_ascii = "".join(
        map(lambda x: chr(x + 33), read_umi_q))
    read_barcode_q_ascii = "".join(
        map(lambda x: chr(x + 33), read_barcode_q))

    # Build alignments and return the parsed results."""
    # Mock a parasail alignemnt result
    # 100 bp prefic sequence
    p_alignment = Mock()
    p_alignment.traceback.query = \
        f"{align_adapter1}{align_barcode}{umi}{polyt}"

    # The probe (reference) would contain Ns at the barcode and UMI
    # positions. The first N will be at position 10 as the
    reference_probe = \
        f"{adapter1_probe_suffix}{'N'* barcode_length}{'N'*umi_length}{polyt}"
    p_alignment.traceback.ref = reference_probe

    (
        adapter1_editdist, barcode_result, umi_result,
        bc_qscores, umi_qscores, bc_min_q
    ) = parse_probe_alignment(
            p_alignment, adapter1_probe_suffix,
            barcode_length, umi_length, prefix_qual)

    assert adapter1_editdist == adapter1_edit_dist
    assert barcode_result == barcode_probe
    assert umi_result == umi
    assert bc_qscores == read_barcode_q_ascii
    assert umi_qscores == read_umi_q_ascii
    assert bc_min_q == min(read_barcode_q)

    # assert adapter1_editdist == 0
    # assert barcode_result == barcode
    # assert umi_result == umi
    # assert bc_qscores == setup.read_barcode_q_ascii
    # assert umi_qscores == setup.read_umi_q_ascii
    # assert bc_min_q == min(setup.read_barcode_q)
    #
    # # no Ns found in reference probe
    # (
    #     adapter1_editdist, barcode_result, umi_result,
    #     bc_qscores, umi_qscores, bc_min_q
    # ) = run(first_n_pos=-1)
    # assert adapter1_editdist == 10
    # assert barcode_result == ""
    # assert umi_result == ""
    # assert bc_min_q == 0
    #
    # # Single substitution in the adapter1 read
    # a1_test = list(adapter1_probe_seq)
    # a1_test[5] = 'A'
    # a1_test = ''.join(a1_test)
    # (
    #     adapter1_editdist, barcode_result, umi_result,
    #     bc_qscores, umi_qscores, bc_min_q
    # ) = run(adapter1_probe_test=a1_test)
    # assert adapter1_editdist == 1
    #
    # # Test a deletion in the read barcode
    # bc_del = list(barcode)
    # bc_del[5] = '-'
    # bc_del = ''.join(bc_del)
    # (
    #     adapter1_editdist, barcode_result, umi_result,
    #     bc_qscores, umi_qscores, bc_min_q
    # ) = run(barcode_test=bc_del)
    # assert adapter1_editdist == 0
    # # The deltion character should have been stripped
    # assert barcode_result == bc_del.replace('-', '')
    # assert umi_result == umi
    # # Check quality scores.
    # assert bc_qscores == setup.read_barcode_q_ascii
    # assert umi_qscores == setup.read_umi_q_ascii
    # assert bc_min_q == min(setup.read_barcode_q)
