"""Test expression matrix construction."""
import os
from pathlib import Path
import tempfile

import pandas as pd
import pytest
from workflow_glue.expression_matrix import main


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
        'CB': ['TTT', 'AAA', 'AAA', 'AAA', 'AAA', 'AAA'],
        'UB': ['CAT', 'AAA', 'CCC', 'TTT', 'AAA', 'GGG']
    })

    expected_raw_result = pd.DataFrame({
        'gene': ['g1', 'g2', 'g3'],
        'TTT': [1, 0, 0],
        'AAA': [1, 1, 2]
    })

    # With a min cells per gene of two, would exclude g1 and g2
    # TODO add more tests

    # Multiply by norm count, divide by total cell counts,
    # np.log1p transform, then divide by log(10) to get back to base 10
    expected_processed_result = pd.DataFrame({
        'gene': ['g1'],
        'TTT': [1.04139],  # np.log1p((1 * 10) / 1) / np.log(10)
        'AAA': [0.54407]   # np.log1p((1 * 10) / 4) / np.log(10)
    })

    return df, expected_raw_result, expected_processed_result


def test_main(tags_df):
    """Test the main function.

    :param tags_df: fixture with input test file and expected result pd.DataFrame.
    :return:
    """
    tags_df, expected_raw_result, expected_processed_result = tags_df
    # df_file = tempfile.NamedTemporaryFile(
    #     suffix='.tsv', mode='w', delete=False).name
    # out_file = tempfile.NamedTemporaryFile(
    #     suffix='.tsv', mode='w', delete=False).name
    with tempfile.TemporaryDirectory() as fh:
        tmp_test_dir = Path(fh)
        tags_file = tmp_test_dir / "tsv1.tsv"
        tags_df.to_csv(tags_file, sep='\t')
        os.chdir(tmp_test_dir)

        class Args:
            read_tags = tags_file
            feature_type = 'gene'
            chunk_size = 1
            output_prefix = 'test_gene'
            min_features = 1
            min_cells = 2
            max_mito = 5
            mito_prefixes = ['MT-']
            norm_count = 10

        main(Args)

        raw_results_df_file = (
                tmp_test_dir / f'{Args.output_prefix}_expression.count.tsv'
        )
        counts_result_df = pd.read_csv(raw_results_df_file, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_raw_result, counts_result_df, check_like=True, check_dtype=False)

        proc_results_df_file = (
                tmp_test_dir / f'{Args.output_prefix}_expression.processed.tsv'
        )
        procs_result_df = pd.read_csv(proc_results_df_file, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_processed_result,
            procs_result_df, check_like=True, check_dtype=False)
