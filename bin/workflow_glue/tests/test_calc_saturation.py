"""Test adapter_scan_vsearch."""
from unittest.mock import Mock

import polars as pl
from workflow_glue.calc_saturation import (
    downsample_dataframe, run_jobs
)


def test_run_jobs(tmp_path):
    """Test_downsample_reads.

    Check for the correct number of downsampled dataframes are returned, each with
    the correct size.
    """
    args = Mock()
    args.read_tags = tmp_path / 'read_tags.tsv'
    args.output = tmp_path / 'output.tsv'
    args.threads = 2
    args.sample = 'test'

    # Create df with 1000 rows of fake data.
    with open(args.read_tags, 'w') as fh:
        fh.write('read_id\tcorrected_barcode\tcorrected_umi\tgene\n')
        row = 'id\tagtcgatcgatcgta\tatcgtacaatct\tYFG'
        for i in range(1000):
            fh.write(f'{row}\n')

    run_jobs(args)
    result = pl.read_csv(source=args.output, separator='\t')

    # Simply check correct number of results are returned
    # and that the downsampled reads are the correct size.
    assert len(result) == 16
    for row in result.iter_rows(named=True):
        assert row['downsamp_reads'] == 1000 * row['downsamp_frac']


def test_downsample_dataframe():
    """Test calc_saturation."""
    header = ['barcode', 'umi', 'gene']

    rows = (
        # Cell 1: 4 reads, 2 umis with two reads each, 2 genes.
        ('AGATAGATAGATAGAT', 'ATAGATAGATAG', 'YFG1'),
        ('AGATAGATAGATAGAT', 'ATAGATAGATAG', 'YFG1'),
        ('AGATAGATAGATAGAT', 'ccccATAGATAG', 'YFG2'),
        ('AGATAGATAGATAGAT', 'ccccATAGATAG', 'YFG2'),

        # Cell 2: 4 reads, 3 umis, 3 genes.
        ('TATATATATATATATA', 'TACTACTACTAC', 'YFG3'),
        ('TATATATATATATATA', 'CACTACTACTCA', 'YFG4'),
        ('TATATATATATATATA', 'CACTACTACTCA', 'YFG4'),
        ('TATATATATATATATA', 'GACGACGACGAC', 'YFG5')
    )

    df = pl.from_records(
        data=rows, schema=header)

    (
        label,
        n_reads,
        reads_per_cell,
        genes_per_cell,
        umis_per_cell,
        umi_saturation
    ) = downsample_dataframe(df, 1.0)

    assert n_reads == 8
    assert reads_per_cell == 4
    assert genes_per_cell == 2.5
    assert umis_per_cell == 2.5

    unique_umis = 5
    n_reads = 8
    assert umi_saturation == 1 - (unique_umis / n_reads)
