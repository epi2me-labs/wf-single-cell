"""Test expression matrix construction."""
import os
from pathlib import Path
import tempfile
from unittest.mock import Mock

import h5py
import numpy as np
import pandas as pd
import pytest
import scipy.sparse
from workflow_glue.expression_matrix import ExpressionMatrix
from workflow_glue.process_matrix import main


@pytest.fixture(params=[False, True], ids=["dense", "sparse"])
def matrix_mode(request):
    """Fixture to provide matrix mode (dense or sparse)."""
    return request.param


@pytest.fixture
def empty_matrix(matrix_mode):
    """Create an empty ExpressionMatrix."""
    matrix = ExpressionMatrix(
        matrix=np.ndarray(shape=(0, 0)),
        features=np.array([]),
        cells=np.array([]),
        sparse=matrix_mode
    )
    return matrix


@pytest.fixture
def small_matrix(matrix_mode):
    """Create a 2x2 ExpressionMatrix."""
    matrix = ExpressionMatrix(
        matrix=np.array([2, 2, 2, 2]).reshape((2, 2)),
        features=np.array(['gene1', 'gene2'], dtype=bytes),
        cells=np.array(['cell1', 'cell2'], dtype=bytes),
        sparse=matrix_mode
    )
    return matrix


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


def test_empty_matrix(matrix_mode):
    """Test instantiating ExpressionMatrix with empty data."""
    em = ExpressionMatrix(
        matrix=np.ndarray(shape=(0, 0)),
        features=np.array([], dtype=bytes),
        cells=np.array([], dtype=bytes),
        sparse=matrix_mode
    )

    assert em.matrix.shape == (0, 0)
    assert em.features.shape == (0,)
    assert em.cells.shape == (0,)


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


def test_aggregate_hdfs(matrix_mode):
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

    em = ExpressionMatrix.aggregate_hdfs(
        (hdf1.name, hdf2.name, hdf3.name), sparse=matrix_mode)
    mat = em.matrix.toarray() if matrix_mode else em.matrix

    np.testing.assert_array_equal(
        em.tcells, np.array(['cell1', 'cell2', 'cell3']))
    np.testing.assert_array_equal(
        em.tfeatures, np.array(['f1', 'f2']))
    np.testing.assert_array_equal(
        mat,  np.array([[5, 2, 2], [5, 4, 1]])
    )


def test_remove_cells_and_features(matrix_mode):
    """Test removing cells and features from ExpressionMatrix."""
    em = ExpressionMatrix(
        matrix=np.array([[0, 5, 0], [1, 1, 1]]),
        features=np.array(['f1', 'f2'], dtype=bytes),
        cells=np.array(['c1', 'c2', 'c3'], dtype=bytes),
        sparse=matrix_mode,
    )

    em.remove_cells(threshold=1)  # no cells should be removed
    em.remove_features(threshold=2)  # only the 5 remains

    np.testing.assert_array_equal(em.tfeatures, np.array(['f2']))
    np.testing.assert_array_equal(em.tcells, np.array(['c1', 'c2', 'c3']))


def test_remove_features_and_cells(matrix_mode):
    """Test removing features and cells from ExpressionMatrix."""
    em = ExpressionMatrix(
        matrix=np.array([[0, 5, 0], [1, 1, 1]]),
        features=np.array(['f1', 'f2'], dtype=bytes),
        cells=np.array(['c1', 'c2', 'c3'], dtype=bytes),
        sparse=matrix_mode,
    )

    em.remove_features(threshold=2)  # only the 5 remains
    em.remove_cells(threshold=1)  # no cells should be removed

    np.testing.assert_array_equal(em.tfeatures, np.array(['f2']))
    np.testing.assert_array_equal(em.tcells, np.array(['c1', 'c2', 'c3']))


def test_remove_cells_then_features_order_sensitive(matrix_mode):
    """Test order-sensitive removal: remove cells first, then features."""
    em = ExpressionMatrix(
        matrix=np.array([
            [1, 0, 0],  # f1
            [1, 1, 0],  # f2
            [0, 1, 1],  # f3
        ]),
        features=np.array(['f1', 'f2', 'f3'], dtype=bytes),
        cells=np.array(['c1', 'c2', 'c3'], dtype=bytes),
        sparse=matrix_mode,
    )

    em.remove_cells(threshold=2)
    em.remove_features(threshold=2)

    np.testing.assert_array_equal(em.tfeatures, np.array(['f2']))
    np.testing.assert_array_equal(em.tcells, np.array(['c1', 'c2']))


def test_remove_features_then_cells_order_sensitive(matrix_mode):
    """Test order-sensitive removal: remove features first, then cells."""
    em = ExpressionMatrix(
        matrix=np.array([
            [1, 0, 0],  # f1
            [1, 1, 0],  # f2
            [0, 1, 1],  # f3
        ]),
        features=np.array(['f1', 'f2', 'f3'], dtype=bytes),
        cells=np.array(['c1', 'c2', 'c3'], dtype=bytes),
        sparse=matrix_mode,
    )

    em.remove_features(threshold=2)
    em.remove_cells(threshold=2)

    np.testing.assert_array_equal(em.tfeatures, np.array(['f2', 'f3']))
    np.testing.assert_array_equal(em.tcells, np.array(['c2']))


def test_normalize_and_log_transform(matrix_mode):
    """Test normalization and log transformation of ExpressionMatrix."""
    em = ExpressionMatrix(
        matrix=np.array([[10, 20], [30, 40]], dtype=float),
        features=np.array(['f1', 'f2'], dtype=bytes),
        cells=np.array(['c1', 'c2'], dtype=bytes),
        sparse=matrix_mode,
    )

    em.normalize(100).log_transform()
    mat = em.matrix.toarray() if matrix_mode else em.matrix

    exp = np.array([
        [1.41497335, 1.53571597],
        [1.88081359, 1.83037478]
    ])
    np.testing.assert_allclose(mat, exp, rtol=1e-5)


def test_mean_and_median_stats(matrix_mode):
    """Test mean and median statistics of ExpressionMatrix."""
    em = ExpressionMatrix(
        matrix=np.array([
            [0, 1, 2],
            [0, 3, 0]
        ]),
        features=np.array(['f1', 'f2'], dtype=bytes),
        cells=np.array(['c1', 'c2', 'c3'], dtype=bytes),
        sparse=matrix_mode,
    )

    # the mean value per cell
    np.testing.assert_array_equal(em.mean_expression, np.array([0, 2, 1]))
    # find total count per cell [0, 4, 2], then take median
    np.testing.assert_array_equal(em.median_counts, np.array([2.]))
    # find non-zero entries per cell [0, 2, 1], then take median
    np.testing.assert_array_equal(em.median_features_per_cell, np.array([1]))


def test_feature_sort(matrix_mode):
    """Test sorting helper attributes are kept up to date."""
    em = ExpressionMatrix(
        matrix=np.array([[0, 2], [3, 4]]),
        features=np.array(['z_gene', 'a_gene'], dtype=bytes),
        cells=np.array(['c2', 'c1'], bytes),
        sparse=matrix_mode,
    )
    assert np.array_equal(em._s_features, [1, 0])
    assert np.array_equal(em._s_cells, [1, 0])

    em.remove_features(threshold=2)
    assert np.array_equal(em._s_cells, [1, 0])
    # f1 is removed, but f2 renumbered
    assert np.array_equal(em._s_features, [0])


def test_cell_sort(matrix_mode):
    """Test sorting helper attributes are kept up to date."""
    em = ExpressionMatrix(
        matrix=np.array([[0, 2], [3, 4]]),
        features=np.array(['z_gene', 'a_gene'], dtype=bytes),
        cells=np.array(['c2', 'c1'], bytes),
        sparse=matrix_mode,
    )
    assert np.array_equal(em._s_features, [1, 0])
    assert np.array_equal(em._s_cells, [1, 0])

    em.remove_cells(threshold=2)
    assert np.array_equal(em._s_cells, [0])
    assert np.array_equal(em._s_features, [1, 0])


def test_main(tags_df):
    """Test the main function of the expression matrix module."""
    tags_df, expected_raw_result, expected_processed_result = tags_df
    with tempfile.TemporaryDirectory() as fh:
        tmp_test_dir = Path(fh)
        os.chdir(tmp_test_dir)
        tags_df.to_csv('tags.tsv', sep='\t')

        args = Mock()
        args.input = ["tags.tsv"]
        args.feature = 'gene'
        args.raw = 'raw.tsv'
        args.per_cell_mito = 'per_cell_mito.tsv'
        args.per_cell_expr = 'per_cell_expr.tsv'
        args.filtered_mex = 'filtered_mex.tsv'
        args.min_features = 1
        args.min_cells = 2
        args.max_mito = 5
        args.mito_prefixes = 'MT-'
        args.norm_count = 10
        args.stats = 'stats.tsv'
        args.processed = 'processed.tsv'
        args.enable_filtering = True
        args.text = True
        args.enable_umap = False
        args.pcn = None

        main(args)

        counts_result_df = pd.read_csv(args.raw, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_raw_result, counts_result_df,
            check_like=True, check_dtype=False)

        procs_result_df = pd.read_csv(args.processed, sep='\t', index_col=None)
        pd.testing.assert_frame_equal(
            expected_processed_result,
            procs_result_df, check_like=True, check_dtype=False)

        mito_results_df = pd.read_csv(
            args.per_cell_mito, sep='\t', index_col=None)
        assert "CB" in mito_results_df.columns
        assert "mito_pct" in mito_results_df.columns


def test_bin_cells_by_coordinates_wrong_cell_names(matrix_mode):
    """Test function raises error on wrong cell names."""
    features = np.array([b"gene1", b"gene2"])
    coords = [(x, y) for x in range(2) for y in range(2)]
    cell_names = [f"b_{x}_{y}".encode() for (x, y) in coords]
    matrix_data = np.array([
        [1, 2, 3, 4],   # gene1
        [5, 6, 7, 8]    # gene2
    ])

    if matrix_mode:
        matrix_data = scipy.sparse.csr_matrix(matrix_data)
    em = ExpressionMatrix(
        matrix=matrix_data, features=features, cells=cell_names, sparse=matrix_mode)

    with pytest.raises(ValueError):
        em.bin_cells_by_coordinates(bin_size=3)


def test_bin_cells_by_coordinates_simple(matrix_mode):
    """Test binning cells by coordinates with a simple example."""
    features = np.array([b"gene1", b"gene2"])
    coords = [(x, y) for x in range(2) for y in range(2)]
    cell_names = [f"s_002um_{x:05d}_{y:05d}-1".encode() for (x, y) in coords]
    matrix_data = np.array([
        [1, 2, 3, 4],   # gene1
        [5, 6, 7, 8]    # gene2
    ])

    if matrix_mode:
        matrix_data = scipy.sparse.csr_matrix(matrix_data)
    em = ExpressionMatrix(
        matrix=matrix_data, features=features, cells=cell_names, sparse=matrix_mode)

    # bin size 2 means all four cells fall into a single bin: (0,0)
    binned = em.bin_cells_by_coordinates(bin_size=2)
    assert binned.matrix.shape == (2, 1)  # 2 genes, 1 bin
    assert binned.sparse == matrix_mode, "Should preserve sparse mode"
    np.testing.assert_array_equal(binned.cells, [b"bin_0_0"])
    expected = matrix_data.sum(axis=1).reshape(2, 1)  # sum all cells
    if matrix_mode:
        binned._matrix = binned.matrix.toarray()
    np.testing.assert_array_equal(binned.matrix, expected)


standard_order = list(range(16))  # simple case, no shuffling
permuted_order = [
    5, 2, 14, 0,
    7, 3, 13, 1,
    10, 6, 15, 4,
    11, 8, 9, 12
]


@pytest.fixture(
    params=[standard_order, permuted_order], ids=["standard", "permuted"])
def permution_order(request):
    """Fixture to provide matrix mode (dense or sparse)."""
    return request.param


def test_bin_cells_by_coordinates_complex(permution_order, matrix_mode):
    """Test binning cells by coordinates with a complex example."""
    # Construct a 4x4 grid of coordinates
    coords = [(x, y) for x in range(4) for y in range(4)]
    cell_names = [f"s_002um_{x:05d}_{y:05d}-1".encode() for (x, y) in coords]

    shuffled_cells = np.array(cell_names)[permution_order]

    # lets create a grid where we expect each resultant to be a sum of
    # 4 original cells with the new cell id. We row-major cell coordinates
    # in the cells name vector, it looks like:
    feature0 = np.array([
        1, 1, 2, 2,
        1, 1, 2, 2,
        3, 3, 4, 4,
        3, 3, 4, 4])[permution_order]  # row major, permuted (maybe)
    matrix = np.array([feature0, feature0 * 10])  # second feature is scaled
    if matrix_mode:
        matrix = scipy.sparse.csr_matrix(matrix)

    em = ExpressionMatrix(
        matrix=matrix,
        features=np.array([b'f1', b'f2']),
        cells=shuffled_cells,
        sparse=matrix_mode
    )
    binned = em.bin_cells_by_coordinates(bin_size=2)

    assert binned.matrix.shape == (2, 4), "Should have 2 features, 4 binned cells"
    assert len(binned.cells) == 4, "Should produce exactly 4 meta-cell bins"
    assert binned.sparse == matrix_mode, "Should preserve sparse mode"

    # the expected matrix is not a function of the permution order,
    # since the binning implicitely reorders the cells by coordinates
    # to always give us back the same
    exp_mat = np.array([
        [4, 8, 12, 16],
        [40, 80, 120, 160]
    ])
    if matrix_mode:
        binned._matrix = binned.matrix.toarray()
    np.testing.assert_array_equal(binned.matrix, exp_mat)
