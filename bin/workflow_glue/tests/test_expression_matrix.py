"""Test expression matrix construction."""
import os
from pathlib import Path
import tempfile

import pandas as pd
import pytest
from workflow_glue.process_matrix import argparser, main


@pytest.fixture()
def tags_df():
    """Make read tag dataframe."""
    # Create a DataFrame with two cells, and genes
    # cell TTT as only one gene (g1) and one UMI (TTT)
    # cell AAA has 3 genes
    #   g1 has a single UMI (TTT)
    #   g2 has a single UMI (AAA)
    #   g3 has two UMIS (CCC, GGG)
    df = pd.DataFrame({
        'gene': ['g1', 'g2', 'g3', 'g1', 'g2', 'g3'],
        'corrected_barcode': ['TTT', 'AAA', 'AAA', 'AAA', 'AAA', 'AAA'],
        'corrected_umi': ['CAT', 'AAA', 'CCC', 'TTT', 'AAA', 'GGG']
    })

    expected_raw_result = pd.DataFrame({
        'gene': ['g1', 'g2', 'g3'],
        'AAA': [1, 1, 2],
        'TTT': [1, 0, 0],
    })

    # With a min cells per gene of two, would exclude g1 and g2
    # TODO add more tests

    # Multiply by norm count, divide by total cell counts,
    # np.log1p transform, then divide by log(10) to get back to base 10
    expected_processed_result = pd.DataFrame({
        'gene': ['g1'],
        'AAA': [1.04139],
        'TTT': [1.04139],  # np.log1p((1 * 10) / 1) / np.log(10)
    })

    return df, expected_raw_result, expected_processed_result


def test_main(tags_df):
    """Test the main function.

    :param tags_df: fixture with input test file and expected result pd.DataFrame.
    :return:
    """
    tags_df, expected_raw_result, expected_processed_result = tags_df
    with tempfile.TemporaryDirectory() as fh:
        tmp_test_dir = Path(fh)
        os.chdir(tmp_test_dir)
        tags_file = "tsv1.tsv"
        tags_df.to_csv(tags_file, sep='\t')

        parser = argparser()
        args = parser.parse_args(
            f"{tags_file} --feature gene --min_features 1 --min_cells 2 "
            "--max_mito 5 --mito_prefixes 'MT-' --norm_count 10 "
            "--enable_filtering --text".split())
        main(args)

        counts_result_df = pd.read_csv(args.raw, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_raw_result, counts_result_df, check_like=True, check_dtype=False)

        procs_result_df = pd.read_csv(args.processed, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_processed_result,
            procs_result_df, check_like=True, check_dtype=False)

        mito_results_df = pd.read_csv(args.per_cell_mito, sep='\t', index_col=None)
        assert "CB" in mito_results_df.columns
        assert "mito_pct" in mito_results_df.columns
