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
    read_tags = test_out_dir / 'sample1.read_summary.tsv'

    assert read_tags.is_file()

    df = pd.read_csv(read_tags, sep='\t')

    # As all reads should be assigned a barcode and UMI, there should be the
    # same number of output rows as in reads in the integration test data (1850).
    assert len(df) == 1850

    # Extract the expected values from the read_id
    df[['true_gene', 'true_transcript', 'true_bc', 'true_umi', 'true_status', '_']] \
        = df['read_id'].str.split('|', expand=True)

    # Check barcode and umis are correctly identified. Allow for 2 incorrect values
    df_barcode_mismatches = df[df.corrected_barcode != df.true_bc]
    assert len(df_barcode_mismatches) < 2

    df_umi_mismatches = df[df.true_umi != df.corrected_umi]
    assert len(df_umi_mismatches) < 2

    # Check gene assignment
    df_gene_matches = df[df.gene == df.true_gene]
    perc_correct = 100 / len(df) * len(df_gene_matches)
    assert perc_correct == 100.0

    # Check transcript assignment
    # We should be getting more than 85% of the transcritps correctly called,
    # especially on this contrived synthetic dataset.
    df_tr_matches = df[df.transcript == df.true_transcript]
    perc_correct = 100 / len(df) * len(df_tr_matches)
    assert perc_correct > 85.0
