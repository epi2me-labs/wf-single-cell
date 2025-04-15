"""Test tag_bam.py"."""
import subprocess as sub

import pandas as pd
import pysam
import pytest
from workflow_glue import tag_bam


def make_bam(tmp_path, bam_entries):
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
    sam_file = tmp_path / 'align.sam'
    with open(sam_file, 'w') as fh_sam:
        fh_sam.write(sam)

    test_bam = tmp_path / 'align.bam'
    sub.check_output(['samtools', 'view', sam_file, '-o', test_bam])
    sub.check_output(['samtools', 'index', test_bam])
    return test_bam


@pytest.mark.parametrize(
    "tags,prim_records,supp_records",
    [
        (  # A single entry from a tag file representing a primary alignment.
            [
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
                {
                    'read_id': 'read2',
                    'CR': 'GCGCGCGCGCGCGCGc',
                    'CB': 'GCGCGCGCGCGCGCCC',
                    'CY': '????????????????',
                    'UR': 'TTTTTTTaTTTT',
                    'UB': 'TTTTTTTTTTTT',
                    'UY': '????????????',
                    'GN': 'YFG2',
                    'TR': 'YFT2',
                    'chr': 'chr9',
                    'start': 1000,
                    'end': 2000
                }
            ],
            # Primary records
            [
                ('read1', 'chr1', 0, "SA:Z:chr2,10000,+,10S100M1S,60,11"),
                ('read2', 'chr9', 0, "")],
            # Supplementary records
            [('read1', 'chr2', 2048, "")]
         )
    ]
)
def test_add_tags(tmp_path, tags, prim_records, supp_records):
    """Check that the output BAMs are tagged correctly."""
    bam_out = tmp_path / 'test_tags.bam'

    tags_to_test = ['CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'GN', 'TR']

    tag_rows = []
    for tag_entry in tags:
        tag_rows.append(pd.DataFrame.from_dict(tag_entry, orient='index').T)
    tags_df = pd.concat(tag_rows, axis=0)
    tags_df.set_index('read_id', drop=True, inplace=True)

    # Get the SA tags; a subset of the primary tags that have a suppl record
    supp_read_ids = [x[0] for x in supp_records]
    sa_tags_df = tags_df.loc[supp_read_ids]

    # Create primary and supplementary tag files in temporary directories
    prim_tags_dir = tmp_path / 'tags'
    prim_tags_dir.mkdir()
    sa_tags_dir = tmp_path / 'sa_tags'
    sa_tags_dir.mkdir()

    # Write the per chr primary tags file
    for chr_, chr_df in tags_df.groupby('chr'):
        chr_df.to_csv(prim_tags_dir / f'{chr_}.tsv', sep='\t')

    sa_tags_file = sa_tags_dir / 'sa_tags.tsv'
    sa_tags_df.to_csv(sa_tags_file, sep='\t')

    test_bam = make_bam(tmp_path, prim_records + supp_records)

    # Run the test
    tag_bam.add_tags(prim_tags_dir, sa_tags_dir, test_bam, bam_out, threads=1)

    # Check that the correct tags have been set on primary and supplementary
    # record.
    primary_tagged = 0
    supp_tagged = 0
    with pysam.AlignmentFile(bam_out, "rb") as bam_result:
        for align in bam_result:
            expected_tags = tags_df.loc[align.query_name]
            for expected_tag, expected_value in expected_tags.items():
                if expected_tag in tags_to_test:
                    assert align.get_tag(expected_tag) == expected_value
            if align.is_supplementary:
                supp_tagged += 1
            else:
                primary_tagged += 1

    assert primary_tagged == len(prim_records)
    assert supp_tagged == len(supp_records)


def test_empty_file(tmp_path):
    """Test giving a header-only tags file, in a tags directory, to tag_bams."""
    input_bam = make_bam(tmp_path, [('read1', 'chr1', 0, '')])
    tags_header = (
        'read_id', 'CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'chr', 'start', 'end', 'gene',
        'transcript')
    tags_df = (
        pd.DataFrame(columns=tags_header)
        .set_index('read_id', drop=True)
        .rename(
            columns={v: k for k, v in tag_bam.BAM_TAGS.items()})
    )
    tmp_test_dir = tmp_path / 'tags'
    tmp_test_dir.mkdir()
    tmp_sa_dir = tmp_path / 'sa_tags'
    tmp_sa_dir.mkdir()
    header_only_file = tmp_test_dir / 'test_tags.tsv'
    header_only_sa_file = tmp_sa_dir / 'test_sa_tags.tsv'
    tags_df.to_csv(header_only_file, sep='\t')
    tags_df.to_csv(header_only_sa_file, sep='\t')
    out_bam = tmp_path / 'out.bam'
    tag_bam.add_tags(tmp_test_dir, tmp_sa_dir, input_bam, out_bam, threads=1)
