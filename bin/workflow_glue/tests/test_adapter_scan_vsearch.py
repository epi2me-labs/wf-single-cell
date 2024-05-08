"""Test adapter_scan_vsearch."""
from pathlib import Path
import tempfile

import pysam
import pytest
from pytest import fixture
from workflow_glue.adapter_scan_vsearch import (
    call_vsearch, complement_trans, create_stranded_reads,
    parse_vsearch, write_adapters_fasta)
from workflow_glue.sc_util import kit_adapters


@fixture
def segment():
    """Random sequence to build a seq by concatenating along with adapters."""
    return (
        "ATTCAGCGCTGGAGACCGAGCGCCCCGCAAAGGGCCTGATCT"
        "ATCGCGCACGGGACTACTCATTGGGACTGCGGCAATAGGGGAGGGGCCTAACAACGTT")


@pytest.mark.parametrize(
    'adapters,expected_results',
    [
        # Non-full length reads
        [[], [['*', 'no_adapters', '*']]],
        [['adapter1_f'], [['adapter1_f', 'single_adapter1', '+']]],
        [['adapter2_r'], [['adapter2_r', 'single_adapter2', '-']]],
        [['adapter2_r', 'adapter1_f'], [['adapter2_r-adapter1_f', 'other', '*']]],
        # Full length reds
        [['adapter1_f', 'adapter2_f'], [['adapter1_f-adapter2_f', 'full_len', '+']]],
        # 3 adapters with one full length segment
        [['adapter2_r', 'adapter1_r', 'adapter1_f'],
         [['adapter2_r-adapter1_r', 'full_len', '-']]],
        # Mutiple subreads in a read
        [['adapter1_f', 'adapter2_f', 'adapter2_r', 'adapter1_r'],
            [
               ['adapter1_f-adapter2_f', 'full_len', '+'],
               ['adapter2_r-adapter1_r', 'full_len', '-'],
        ]],
    ]
)
def test_call_vsearch(adapters, expected_results, segment):
    """Test call_vsearch running and parsing.

    This is the main function of the script that calls a bunch of other functions.
    """
    id_ = 'read_1'
    kits = ['3prime', '5prime', 'multiome']

    for kit in kits:
        # Build the read
        adapter_seqs = []
        for a in adapters:
            # Get name and orientation from eg: adapter2_r
            adapter_name, ori = a.split('_')
            adap = kit_adapters[kit][adapter_name]
            if ori == 'r':
                adap = adap[::-1].translate(complement_trans)
            adapter_seqs.append(adap)

        seq = segment.join(adapter_seqs) + segment

        fastq = (
            f"@{id_}\n"
            f"{seq}\n"
            "+\n"
            f"{'<' * len(seq)}")

        fastq_file = tempfile.NamedTemporaryFile(suffix='.fq')
        with open(fastq_file.name, 'w') as fh:
            fh.write(fastq)

        adapter_fasta = 'adapter_seqs.fasta'
        write_adapters_fasta(
            kit_adapters[kit]['adapter1'], kit_adapters[kit]['adapter2'],
            adapter_fasta)

        vsearch_results = tempfile.NamedTemporaryFile(suffix='.fq')
        call_vsearch(
            Path(fastq_file.name), Path(vsearch_results.name), 0.7, adapter_fasta, 4)
        parsed_results = parse_vsearch(vsearch_results.name)

        # Each result can contain 0 or more subreads -
        #  segments with consecutive pairs of compatible adapters.
        parsed_results = iter(parsed_results[id_])
        for exp_result in expected_results:
            subread_result = next(parsed_results)

            assert subread_result['adapter_config'] == exp_result[0]
            assert subread_result['lab'] == exp_result[1]
            assert subread_result['orig_strand'] == exp_result[2]


def test_write_stranded_fastq():
    """Test that the  correct stranded and trimmed fastq files are being written."""
    # Build a dummy fastq file containing a single read with two subreads.

    # This is a 3prime read
    seq = 't' * 10 + 'A' * 100 + 't' * 10 + 'G' * 200
    fastq = (
        f"@read_1\n"
        f"{seq}\n"
        "+\n"
        f"{'<' * len(seq)}")

    # This config defines one read containing two subreads.
    config = {
        'read_1': [
            {
                'readlen': 100, 'read_id': 'read_1_0', 'start': 10,
                'end': 110,
                'fl': True, 'stranded': True, 'orig_strand': '+',
                'orig_adapter_config':
                    'adapter1_f-adapter2_f-adapter2_r-adapter1_r',
                'adapter_config': 'adapter1_f-adapter2_f',
                'lab': 'full_len'},
            {
                'readlen': 200, 'read_id': 'read_1_1', 'start': 120,
                'end': 320,
                'fl': True, 'stranded': True, 'orig_strand': '-',
                'orig_adapter_config':
                    'adapter1_f-adapter2_f-adapter2_r-adapter1_r',
                'adapter_config': 'adapter1_f-adapter2_f',
                'lab': 'full_len'}
            ]
    }

    temp_fq = tempfile.NamedTemporaryFile(suffix='.fq')

    temp_fq_out = tempfile.NamedTemporaryFile(mode='wt', suffix='.fq')

    with open(temp_fq.name, 'w') as fh:
        fh.write(fastq)

    data = create_stranded_reads(temp_fq.name, config, '3prime', fl_only=True)
    for read in data:
        temp_fq_out.write(read)
    temp_fq_out.flush()

    results = []
    with pysam.FastxFile(temp_fq_out.name) as fh_res:
        for entry in fh_res:
            results.append(entry)
    assert len(results) == 2
    assert len(results[0].sequence) == 100
    # Subread 0 should have been reverse complemented and be all 'T'
    assert set(results[0].sequence) == {'T'}
    assert len(results[1].sequence) == 200
    # Subread 1 should remain unchanged
    assert set(results[1].sequence) == {'G'}
