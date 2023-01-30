"""Integration testing of the whole workflow using synthetic data."""

from pathlib import Path

import pandas as pd
from pytest import fixture


@fixture
def wf_out_dir(request):
    """Set workflow directory."""
    return request.config.getoption('--wf_out_dir')


@fixture
def sample_id(request):
    """Set sample ID."""
    return request.config.getoption('--sample_id')


def test_workflow(wf_out_dir, sample_id):
    """Test the whole Nextflow workflow."""
    out_dir = Path(wf_out_dir)
    test_out_dir = out_dir / sample_id
    read_tags = test_out_dir / 'sample1_read_tags.tsv'

    assert read_tags.is_file()

    df = pd.read_csv(read_tags, sep='\t')

    # Extract the expected values from the read_id
    df[['true_gene', 'true_transcript', 'true_bc', 'true_umi', 'true_status', '_']] \
        = df['read_id'].str.split('|', expand=True)

    # Check barcode and umis are correctly identified
    df_barcode_mismatches = df[df.barcode != df.true_bc]
    assert (len(df_barcode_mismatches) < 2)

    df_umi_mismatches = df[df.true_umi != df.umi]
    assert (len(df_umi_mismatches) < 2)
