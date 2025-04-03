"""Test tag_bam.py"."""
from pathlib import Path
import subprocess as sub
import tempfile

import pandas as pd
import pysam
import pytest
from workflow_glue import tag_bam


def make_bam(bam_entries):
    """Make a BAM file."""
    read = (
        "ATTCAGCGCTGGAGACCGAGCGCCCCGCAAAGGGCCTGATCT"
        "ATCGCGCACGGGACTACTCATTGGGACTGCGGCAATAGGGGAGGGGCCTAACAACGTT")
    chrs = set([x[1] for x in bam_entries])

    # Create the BAM file to be tagged
    header = '\n'.join([f'@SQ	SN:{chr_}	LN:10000000' for chr_ in chrs])

    entries = [f'{header}']

    for records in bam_entries:
        # Make a sam file containing the read and a quality qscore of 60.
        id_, chr_, flag, sa_tag = records
        entries.append(
            f"{id_}\t{flag}\t{chr_}\t1\t60\t{len(read)}M\t*\t0\t0\t"
            f"{read}\t{'?' * len(read)}\t{sa_tag}"
        )
    sam = '\n'.join(entries)
    # Write out a test BAM
    with tempfile.NamedTemporaryFile(mode='w', suffix='.sam', delete=False) as fh_sam:
        fh_sam.write(sam)
        sam_file = fh_sam.name

    test_bam = tempfile.NamedTemporaryFile('w', delete=False, suffix='.bam').name
    sub.check_output(['samtools', 'view', sam_file, '-o', test_bam])
    sub.check_output(['samtools', 'index', test_bam])
    return test_bam


@pytest.mark.parametrize(
    "tags,bam_entries",
    [
        (  # A single entry from a tag file representing a primary alignment.
            {
                'read_id': 'read1',
                'CR': 'AAAAAAAAAAAAAgAA',
                'CB': 'AAAAAAAAAAAAAaAA',
                'CY': '????????????????',
                'UR': 'GGGGGtGGGGGG',
                'UB': 'GGGGGGGGGGGG',
                'UY': '????????????',
                'GN': 'YFG',
                'TR': 'YFT',
                'chr': 'chr1',
                'start': 1000,
                'end': 2000

            },
            # Two alignment records one primary and one supplementary.
            (
                ('read1', 'chr1', 0, "SA:Z:chr2,10000,+,10S100M1S,60,11"),
                ('read1', 'chr2', 2048, "")
            )
         )
    ]
)
def test_add_tags(tags, bam_entries):
    """Check that the output - bams and tag file - are correct."""
    bam_out = tempfile.NamedTemporaryFile('w', suffix='.bam', delete=False).name

    tags_to_test = ['CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'GN', 'TR']
    # Create the tags file
    tags_header = tags.keys()
    tags_rows = [tags.values()]

    tags_df = (
        pd.DataFrame(tags_rows, columns=tags_header)
        .set_index('read_id', drop=True)
        .rename(  # Convert to long tag form - those expected in the tags file
            columns={v: k for k, v in tag_bam.BAM_TAGS.items()})
    )

    tags_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv', delete=False)
    tags_df.to_csv(tags_file.name, sep='\t')

    test_bam = make_bam(bam_entries)

    # Run the test
    tag_bam.add_tags(
        Path(tags_file.name), test_bam, bam_out, threads=1)

    # Check that the correct tags have been set on primary and supplementary
    # record.
    tagged = 0
    with pysam.AlignmentFile(bam_out, "rb") as bam_result:
        for align in bam_result:
            for expected_tag, expected_value in tags.items():
                if expected_tag in tags_to_test:
                    assert align.get_tag(expected_tag) == expected_value
            tagged += 1
    assert tagged == 2


def test_empty_file():
    """Test giving a header-only tags file, in a tags directory, to tag_bams."""
    input_bam = make_bam([('read1', 'chr1', 0, '')])
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
        out_bam_p = tempfile.NamedTemporaryFile('w', suffix='.bam', delete=False).name
        tag_bam.add_tags(tmp_test_dir, input_bam, out_bam_p, threads=1)
