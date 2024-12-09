"""Test tag_bam.py"."""
from pathlib import Path
import subprocess as sub
import tempfile

import pandas as pd
import pysam
import pytest
from workflow_glue import tag_bam


@pytest.fixture
def read():
    """Get a randomly-generated gene sequence."""
    return (
        "ATTCAGCGCTGGAGACCGAGCGCCCCGCAAAGGGCCTGATCT"
        "ATCGCGCACGGGACTACTCATTGGGACTGCGGCAATAGGGGAGGGGCCTAACAACGTT")


@pytest.fixture
def input_bam(read):
    """Create a dummy bam file.

    Each entry has the same the same random sequence but a unique read_id.
    """
    header = """@SQ	SN:chr17	LN:10000000"""
    read_ids = ['read1', 'read2']
    entries = [f'{header}']

    for read_id in read_ids:
        # Make a sam file containing the read and a quality qscore of 60.
        entries.append((
            f"{read_id}\t0\tchr17\t1\t60\t{len(read)}M\t*\t0\t0\t"
            f"{read}\t{'?' * len(read)}"
        ))
    sam = '\n'.join(entries)
    # Write out a test BAM
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as fh_sam:
        fh_sam.write(sam)
        sam_file = fh_sam.name

    bam = tempfile.NamedTemporaryFile('w', delete=False, suffix='.bam').name
    sub.check_output(['samtools', 'view', sam_file, '-o', bam])
    sub.check_output(['samtools', 'index', bam])

    return bam


@pytest.fixture
def tags_file():
    """Create a tags TSV file."""
    tags_header = (
        'read_id', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'chr', 'start', 'end', 'gene',
        'transcript')
    tags_rows = (
        ('read1', 'AAAAAAAAAAAAAgAA', 'AAAAAAAAAAAAAaAA', '????????????????',
         'GGGGGtGGGGGG', 'GGGGGGGGGGGG', '????????????', 'chr17', 1000, 2000,
         'YFG', 'YFT'),
        ('read2', 'TcTTTTTTTTTTTTTT', 'TTTTTTTTTTTTTTTT', '****************',
         'CCCCCCaCCCCC', 'CCCCCCCCCCCC', '????????????', 'chr17', 4000, 5000,
         'YFG1', 'YFT2')

    )
    tags_df = (
        pd.DataFrame(tags_rows, columns=tags_header)
        .set_index('read_id', drop=True)
        .rename(
            columns={v: k for k, v in tag_bam.BAM_TAGS.items()})
    )
    tags = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
    tags_df.to_csv(tags.name, sep='\t')

    return tags.name


def test_add_tags(tags_file, input_bam):
    """Check that the output - bams and tag file - are correct."""
    # create_inputs
    out_bam = tempfile.NamedTemporaryFile('w', suffix='.bam').name
    tag_bam.add_tags(Path(tags_file), input_bam, out_bam, 1)

    # Check that the correct tags have been set
    with pysam.AlignmentFile(out_bam, "rb") as bam_result:
        for align in bam_result:
            read_id = align.query_name
            # Check that some of the tags have been set correctly
            if read_id == 'read1':
                assert align.get_tag('CR') == 'AAAAAAAAAAAAAgAA'
                assert align.get_tag('CB') == 'AAAAAAAAAAAAAaAA'
                assert align.get_tag('CY') == '????????????????'
                assert align.get_tag('UR') == 'GGGGGtGGGGGG'
                assert align.get_tag('UB') == 'GGGGGGGGGGGG'
                assert align.get_tag('UY') == '????????????'


def test_empty_file(input_bam):
    """Test giving a header-only tags file, in a tags directory, to tag_bams."""
    tags_header = (
        'read_id', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'chr', 'start', 'end', 'gene',
        'transcript')
    tags_df = (
        pd.DataFrame(columns=tags_header)
        .set_index('read_id', drop=True)
        .rename(
            columns={v: k for k, v in tag_bam.BAM_TAGS.items()})
    )
    with tempfile.TemporaryDirectory() as fh:
        tmp_test_dir = Path(fh)
        header_only_file = tmp_test_dir / 'test_tags.tsv'
        tags_df.to_csv(header_only_file, sep='\t')
        out_bam = tempfile.NamedTemporaryFile('w', suffix='.bam', delete=False).name
        tag_bam.add_tags(tmp_test_dir, input_bam, out_bam, 1)
