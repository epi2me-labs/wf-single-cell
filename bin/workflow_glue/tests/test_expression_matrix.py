"""Test expression matrix construction."""
import os
from pathlib import Path
import tempfile

import h5py
import numpy as np
import pandas as pd
import pytest
from workflow_glue.expression_matrix import ExpressionMatrix
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


def test_empty_matrix():
    """Test instantiating ExpressionMatrix with empty data."""
    em = ExpressionMatrix(
        matrix=np.ndarray(shape=(0, 2)),
        features=np.array([], dtype=bytes),
        cells=np.array([], dtype=bytes)
    )

    assert em.matrix.shape == (0, 2)
    assert em.features.shape == (0,)
    assert em.cells.shape == (0,)


@pytest.fixture
def empty_matrix():
    """Create an empty ExpressionMatrix."""
    matrix = ExpressionMatrix(
        matrix=np.ndarray(shape=(0, 2)),
        features=np.array([]),
        cells=np.array([])
    )
    return matrix


@pytest.fixture
def small_matrix():
    """Create a 2x2 ExpressionMatrix."""
    matrix = ExpressionMatrix(
        matrix=np.array([2, 2, 2, 2]).reshape((2, 2)),
        features=np.array(['gene1', 'gene2'], dtype=bytes),
        cells=np.array(['cell1', 'cell2'], dtype=bytes)
    )
    return matrix


def test_normalize_empty_matrix(empty_matrix):
    """Test normalizing with empty ExpressionMatrix."""
    empty_matrix.normalize(10000)


def test_remove_cells_empty_matrix(empty_matrix):
    """Test cell filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="Matrix is zero-sized on entry to `remove_cells`."
    ):
        empty_matrix.remove_cells(0)


def test_remove_cells_becomes_empty(small_matrix):
    """Test cell filtering with small matrix made empty by filtering."""
    with pytest.raises(
        ValueError,
        match="All cells would be removed, try altering filter thresholds."
    ):
        small_matrix.remove_cells(3)


def test_remove_features_empty_matrix(empty_matrix):
    """Test feature filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="Matrix is zero-sized on entry to `remove_features`."
    ):
        empty_matrix.remove_features(0)


def test_remove_features_becomes_empty(small_matrix):
    """Test fature filtering with small matrix made empty by filtering."""
    with pytest.raises(
        ValueError,
        match="All features would be removed, try altering filter thresholds."
    ):
        small_matrix.remove_features(3)


def test_remove_cells_and_features_empty_matrix(empty_matrix):
    """Test cell and feature filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="Matrix is zero-sized on entry to `remove_cells_and_features`."
    ):
        empty_matrix.remove_cells_and_features(0, 0)


def test_remove_cells_and_features_becomes_empty(small_matrix):
    """Test cell and feature filtering made empty by filtering."""
    with pytest.raises(
        ValueError,
        match="All features would be removed, try altering filter thresholds."
    ):
        small_matrix.remove_cells_and_features(3, 3)


def test_remove_skewed_cells_empty_matrix(empty_matrix):
    """Test skewed cell filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="Matrix is zero-sized on entry to `remove_skewed_cells`."
    ):
        empty_matrix.remove_skewed_cells(0, ['gene'])


def test_remove_skewed_cells_becomes_empty(small_matrix):
    """Test skewed cell filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="All cells would be removed, try altering filter thresholds."
    ):
        small_matrix.remove_skewed_cells(0.05, ['gene'])


def test_remove_unknown_empty_matrix(empty_matrix):
    """Test unknown feature filtering with empty ExpressionMatrix."""
    with pytest.raises(
        ValueError,
        match="Matrix is zero-sized on entry to `remove_unknown`."
    ):
        empty_matrix.remove_unknown('-')


def test_remove_unknown_becomes_empty():
    """Test unknown filtering made empty by filtering."""
    matrix = ExpressionMatrix(
        matrix=np.array([2, 2]).reshape((1, 2)),
        features=np.array(['-'], dtype=bytes),
        cells=np.array(['cell1', 'cell2'], dtype=bytes)
    )
    with pytest.raises(
        ValueError,
        match="All features would be removed, try altering filter thresholds."
    ):
        matrix.remove_unknown('-')


def test_aggregate_hdfs():
    """Test the creation of an ExpressionMatrix from multiple HDF inputs."""
    hdf1 = tempfile.NamedTemporaryFile(suffix='.hdf5', mode='w')
    with h5py.File(hdf1.name, 'w') as fh1:
        fh1['cells'] = ['cell1', 'cell2']
        fh1['features'] = ['f1', 'f2']
        fh1['matrix'] = np.array([
            [1, 2], [3, 4]
        ]).reshape((2, 2))

    hdf2 = tempfile.NamedTemporaryFile(suffix='.hdf5', mode='w')
    with h5py.File(hdf2.name, 'w') as fh2:
        fh2['cells'] = ['cell3', 'cell1']
        fh2['features'] = ['f2', 'f1']
        fh2['matrix'] = np.array([
            [1, 2], [2, 4]
        ]).reshape((2, 2))

    hdf3 = tempfile.NamedTemporaryFile(suffix='.hdf5', mode='w')
    with h5py.File(hdf3.name, 'w') as fh3:
        fh3['cells'] = []
        fh3['features'] = []
        fh3['matrix'] = np.ndarray(shape=(0, 0))

    em = ExpressionMatrix.aggregate_hdfs((hdf1.name, hdf2.name, hdf3.name))

    np.testing.assert_array_equal(em.tcells, np.array(['cell1', 'cell2', 'cell3']))
    np.testing.assert_array_equal(em.tfeatures, np.array(['f1', 'f2']))
    np.testing.assert_array_equal(
        em.matrix,  np.array([[5, 2, 2], [5, 4, 1]])
    )


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
