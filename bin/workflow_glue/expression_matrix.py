"""Expression counts matrix construction."""
from collections import defaultdict
import copy
import gc
import gzip
import os
import re

import h5py
import numpy as np
import pandas as pd
import polars as pl
import scipy.io
import scipy.sparse

from .util import get_named_logger  # noqa: ABS101

VISIUM_HD_REGEX = re.compile(rb"^s_002um_\d{5}_\d{5}-1$")


class ExpressionMatrix:
    """Representation of expression matrices."""

    def __init__(
            self, matrix=None, features=None, cells=None,
            fname=None, dtype=int, sparse=False, saturation=None):
        """Create a matrix.

        Do not use the constructor directly.
        """
        self.logger = get_named_logger("ExpressionMatrix")
        self._matrix = matrix
        self._fh = None
        self._features = features
        self._cells = cells
        self._dtype = dtype
        self._sparse = sparse
        self._s_features = None
        self._s_cells = None
        self._saturation = saturation

        if fname is not None:
            self._fh = h5py.File(fname, 'r')
        else:
            if any(x is None for x in (features, cells)):
                raise ValueError(
                    "Either fname or features, and cells must be provided.")
            if matrix is None:
                # make an empty matrix
                if not sparse:
                    matrix = np.zeros((len(features), len(cells)), dtype=dtype)
                else:
                    matrix = scipy.sparse.csr_matrix(
                        (len(features), len(cells)), dtype=dtype)
            else:
                if not scipy.sparse.issparse(matrix) and sparse:
                    matrix = scipy.sparse.csr_matrix(
                        matrix, shape=(len(features), len(cells)), dtype=dtype)

            self._matrix = matrix
            self._features = features
            self._cells = cells
            self._s_features = np.argsort(features)
            self._s_cells = np.argsort(cells)
            self._saturation = saturation
            if self._matrix.shape[0] != len(features):
                raise ValueError("Matrix row count does not match number of features.")
            if self._matrix.shape[1] != len(cells):
                raise ValueError("Matrix column count does not match number of cells.")

    def __exit__(self):
        """Cleanup."""
        if self._fh is not None:
            self._fh.close()

    @classmethod
    def from_hdf(cls, name, dtype=int, sparse=False):
        """Load a matrix from HDF file."""
        ma = cls(fname=name, dtype=dtype, sparse=sparse)
        return ma

    @classmethod
    def from_tags(cls, df, feature='gene'):
        """Create a matrix from a tags TSV file, or dataframe."""
        if not isinstance(df, pd.DataFrame):
            df = pd.read_csv(
                df,
                index_col=None,
                usecols=['corrected_barcode', 'corrected_umi', feature],
                sep='\t',
                dtype='category')

        # calculation saturation before dropping unassigned
        saturation = Saturation.from_tags(df)

        df.drop(df.index[df[feature] == '-'], inplace=True)
        df = (
            df.groupby([feature, 'corrected_barcode'], observed=True)
            .nunique()['corrected_umi']
            .reset_index()  # Get feature back to a column
            .pivot(index=feature, columns='corrected_barcode', values='corrected_umi')
        )
        df.fillna(0, inplace=True)

        em = cls(
            df.to_numpy(dtype=int),
            df.index.to_numpy(dtype=bytes),
            df.columns.to_numpy(dtype=bytes),
            dtype=int,
            saturation=saturation
        )
        return em

    @classmethod
    def aggregate_hdfs(cls, fnames, dtype=int, sparse=None):
        """Aggregate a set of matrices stored in HDF."""
        logger = get_named_logger("ExpressionMatrix")

        # bootstrap saturation, and get union of all barcodes and features
        saturation = Saturation.zero()
        features = set()
        cells = set()
        for fname in fnames:
            ma = cls.from_hdf(fname, dtype=dtype)
            features.update(ma.features)
            cells.update(ma.cells)

        # use sparse if specified, otherwise force use if matrix >8 Gb
        estimated_bytes = len(features) * len(cells) * np.dtype(dtype).itemsize
        if sparse is not None:
            use_sparse = sparse
        else:
            use_sparse = estimated_bytes > 8 * (1024 ** 3)
        logger.info(
            f"Aggregating {len(fnames)} matrices, "
            f"estimated size: {estimated_bytes / (1024 ** 3):.2f} Gb, "
            f"sparse: {use_sparse}.")

        # sort by names, just to be nice, doesn't guarantee matrices can be compared
        full_matrix = cls(
            features=np.array(sorted(features), dtype=bytes),
            cells=np.array(sorted(cells), dtype=bytes),
            dtype=dtype, sparse=use_sparse, saturation=saturation)

        for i, fname in enumerate(fnames):
            ma = cls.from_hdf(
                fname, dtype=dtype, sparse=use_sparse)
            full_matrix += ma

        if use_sparse:
            full_matrix._matrix = full_matrix.matrix.tocsr()
        return full_matrix

    @classmethod
    def aggregate_tags(cls, fnames, feature='gene', dtype=int):
        """Aggregate a set of tags files."""
        if len(fnames) == 1:
            return cls.from_tags(fnames[0])
        for fname in fnames:
            ma = ExpressionMatrix.from_tags(fname, feature=feature)
            ma.to_hdf(f"{fname}.hdf")
        hdfs = [f"{fname}.hdf" for fname in fnames]
        return cls.aggregate_hdfs(hdfs, dtype=dtype)

    def to_hdf(self, fname):
        """Save matrix to HDF."""
        with h5py.File(fname, 'w') as fh:
            fh['cells'] = self.cells
            fh['features'] = self.features
            fh['matrix'] = self.matrix
            fh["saturation"] = self.saturation._data

    def to_mex(
            self, fname, feature_type="Gene Expression", feature_ids=None, dtype=None):
        """Export to MEX folder."""
        os.mkdir(fname)
        # barcodes, write bytes directly
        with gzip.open(os.path.join(fname, "barcodes.tsv.gz"), 'wb') as fh:
            for bc in self.cells:
                fh.write(bc)
                fh.write("-1\n".encode())
        # features, convert to str, write as text
        features = self.tfeatures
        if feature_ids is None:
            feature_ids = {x: f"unknown_{i:05d}" for i, x in enumerate(features)}
        with gzip.open(os.path.join(fname, "features.tsv.gz"), 'wt') as fh:
            for feat in features:
                fh.write(f"{feature_ids[feat]}\t{feat}\t{feature_type}\n")
        # matrix as bytes
        with gzip.open(os.path.join(fname, "matrix.mtx.gz"), 'wb') as fh:
            coo = scipy.sparse.coo_matrix(self.matrix, dtype=dtype)
            scipy.io.mmwrite(
                fh, coo,
                comment=(
                    'metadata_json: {'
                    '"software_version": "ont-single-cell",'
                    '"format_version": 2}'))

    # TODO: could store feature type in class
    def to_tsv(self, fname, index_name):
        """Write matrix to tab-delimiter file."""
        self.write_matrix(
            fname, self.matrix, self.tfeatures, self.tcells, index_name=index_name)

    @property
    def features(self):
        """An array of feature names."""
        if self._features is not None:
            return self._features

        if self._fh is not None:
            self._features = self._fh['features'][()].astype(bytes)
            self._s_features = np.argsort(self._features)
        else:
            raise RuntimeError("Features not initialized.")
        return self._features

    @property
    def cells(self):
        """An array of cell names."""
        if self._cells is not None:
            return self._cells
        if self._fh is not None:
            self._cells = self._fh['cells'][()].astype(bytes)
            self._s_cells = np.argsort(self._cells)
        else:
            raise RuntimeError("Cells not initialized.")
        return self._cells

    @property
    def saturation(self):
        """Get saturation data."""
        if self._saturation is not None:
            return self._saturation
        if self._fh is not None:
            self._saturation = Saturation.from_arr(
                self._fh['saturation'][()])
        else:
            raise RuntimeError("saturation not initialized")
        return self._saturation

    @property
    def matrix(self):
        """Expression matrix array."""
        if self._matrix is not None:
            return self._matrix

        if self._fh is not None:
            matrix = self._fh['matrix'].astype(self._dtype)[()]
            if self._sparse:
                matrix = scipy.sparse.csr_matrix(
                    matrix, shape=(len(self.features), len(self.cells)))
            self._matrix = matrix
        else:
            raise RuntimeError("Matrix not initialized.")
        return self._matrix

    @property
    def sparse(self):
        """Return True if this is a sparse matrix."""
        return scipy.sparse.issparse(self._matrix)

    @property
    def is_visium_hd(self):
        """Return True if this is a Visium HD matrix."""
        if VISIUM_HD_REGEX.match(self._cells[0]):
            return True
        else:
            return False

    @property
    def tcells(self):
        """A copy of the cell names as a text array."""
        return np.array(
            [x.decode() for x in self.cells], dtype=str)

    @property
    def tfeatures(self):
        """A copy of the feature names as a text array."""
        return np.array(
            [x.decode() for x in self.features], dtype=str)

    @property
    def mean_expression(self):
        """The mean expression per cell.

        This is a vector of the mean count for each cell,
        averaged across all features.
        """
        if self.sparse:
            sums = np.asarray(self.matrix.sum(axis=0)).flatten()
            mean = sums / self.matrix.shape[0]
        else:
            mean = self.matrix.mean(axis=0)
        return mean

    @property
    def median_counts(self):
        """Median total UMI counts per cell.

        This is the median of the total expression per cell,
        i.e. the sum across features.
        """
        if self.sparse:
            counts = np.asarray(self.matrix.sum(axis=0)).flatten()
        else:
            counts = self.matrix.sum(axis=0)
        return np.median(counts)

    @property
    def median_features_per_cell(self):
        """Median features detected per cell.

        This is the median of the num. of unique (non-zero) features per cell.
        """
        if self.sparse:
            # see scipy.sparse.csr_matrix.count_nonzero.html
            # matrix should be csr already, so tocsr() is no-op but included
            # here for safety
            nonzero = np.bincount(
                self.matrix.tocsr().indices, minlength=self.matrix.shape[1])
        else:
            nonzero = np.count_nonzero(self.matrix, axis=0)
        return np.median(nonzero)

    @property
    def saturation_statistics(self):
        """Get feature and UMI statistics at varying sampling fractions."""
        mat = copy.copy(self)  # maybe will need to do in place for memory?
        mat._matrix = mat.matrix.astype(int)
        arr = np.zeros(
            len(Saturation.FRACTIONS),
            dtype=[
                ("frac", float), ("n_umis", int),
                ("median_feats", int), ("median_umis", int)])
        last_frac = None
        for i, frac in enumerate(Saturation.FRACTIONS):
            p = frac if last_frac is None else frac / last_frac
            mat.downsample(p)
            arr["frac"][i] = frac
            arr["n_umis"][i] = np.sum(mat.matrix)
            arr["median_feats"][i] = mat.median_features_per_cell
            arr["median_umis"][i] = mat.median_counts
            last_frac = frac
        return pd.DataFrame(arr)

    def normalize(self, norm_count):
        """Normalize total cell weight to fixed cell count."""
        total_per_cell = np.asarray(self.matrix.sum(axis=0)).flatten()
        scaling = norm_count / total_per_cell
        if self.sparse:
            self._matrix = self.matrix.multiply(scaling)
        else:
            self._matrix = np.multiply(self.matrix, scaling, dtype=float)
        return self

    def remove_cells(self, threshold):
        """Remove cells with few features present."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_cells`.")

        if self.sparse:
            n_features = self.matrix.getnnz(axis=0)
        else:
            n_features = np.count_nonzero(self.matrix, axis=0)

        mask = n_features >= threshold
        if np.count_nonzero(mask) == 0:
            raise ValueError(
                "All cells would be removed, try altering filter thresholds.")
        self._remove_elements(cell_mask=mask)
        return self

    def remove_features(self, threshold):
        """Remove features that are present in few cells."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_features`.")

        if self.sparse:
            n_cells = self.matrix.getnnz(axis=1)
        else:
            n_cells = np.count_nonzero(self.matrix, axis=1)

        mask = n_cells >= threshold
        if np.count_nonzero(mask) == 0:
            raise ValueError(
                "All features would be removed, try altering filter thresholds.")
        self._remove_elements(feat_mask=mask)
        return self

    def remove_cells_and_features(self, cell_thresh, feature_thresh):
        """Remove cells and features simultaneously.

        The feature and cell removal methods do not commute, this
        method provide a more uniform handling of cell and feature
        filtering.
        """
        if self.matrix.size == 0:
            raise ValueError(
                "Matrix is zero-sized on entry to `remove_cells_and_features`.")

        if self.sparse:
            n_features = self.matrix.getnnz(axis=0)
            n_cells = self.matrix.getnnz(axis=1)
        else:
            n_features = np.count_nonzero(self.matrix, axis=0)
            n_cells = np.count_nonzero(self.matrix, axis=1)

        cell_mask = n_features >= cell_thresh
        feat_mask = n_cells >= feature_thresh
        self._remove_elements(feat_mask=feat_mask, cell_mask=cell_mask)

        return self

    def remove_skewed_cells(self, threshold, prefixes, fname=None, label=None):
        """Remove cells with overabundance of feature class."""
        if self.matrix.size == 0:
            raise ValueError(
                "Matrix is zero-sized on entry to `remove_skewed_cells`.")
        sel = self.find_features(prefixes)
        sel_total = self.matrix[sel].sum(axis=0, dtype=float)
        total = self.matrix.sum(axis=0, dtype=float)
        sel_total /= total
        mask = sel_total <= threshold
        if fname is not None:
            # providing prefixes is a small list, this won't be too much data
            self.logger.info("Writing skewed cell percentages to %s", fname)
            self.write_matrix(
                fname, 100 * sel_total, self.tcells, [label], index_name="CB")
        self._remove_elements(cell_mask=mask)
        return self

    def remove_unknown(self, key="-"):
        """Remove the "unknown" feature."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_unknown`.")
        mask = self.features == key.encode()
        self._remove_elements(feat_mask=~mask)

    def find_features(self, prefixes, inverse=False):
        """Find features by name prefix."""
        sel = np.array([], dtype=int)
        feat = self.tfeatures
        for pre in prefixes:
            sel = np.union1d(sel, np.argwhere(np.char.find(feat, pre) > -1))
        if inverse:
            sel = np.setdiff1d(np.arange(len(self.features)), sel)
        return sel

    def log_transform(self):
        """Transform expression matrix to log scale."""
        m = self._matrix.data if self.sparse else self._matrix
        np.log1p(m, out=m)
        m /= np.log(10)
        return self

    def downsample(self, p):
        """Binomial down sampling of the counts matrix."""
        if self.sparse:
            self.matrix.data = np.random.binomial(n=self.matrix.data, p=p)
            self.matrix.eliminate_zeros()
        else:
            self._matrix = np.random.binomial(n=self.matrix, p=p)
        return self

    @staticmethod
    def write_matrix(fname, matrix, index, col_names, index_name="feature"):
        """Write a matrix (or vector) to TSV."""
        # Ensure matrix is dense! This might be bad
        if scipy.sparse.issparse(matrix):
            matrix = matrix.toarray()
        matrix = np.atleast_2d(matrix)

        # transpose if the shape doesn't match expected dimensions
        if matrix.shape[0] != len(index):
            if matrix.shape[1] == len(index) and matrix.shape[0] == len(col_names):
                matrix = matrix.T
            else:
                raise ValueError(
                    f"Shape of matrix {matrix.shape} is incompatible with "
                    f"index (len {len(index)}) and columns (len {len(col_names)})."
                )
        df = pd.DataFrame(matrix, copy=False, columns=col_names, index=index)
        df.index.rename(index_name, inplace=True)
        df = pl.from_pandas(df, include_index=True)
        df.write_csv(fname, separator="\t")
        del df
        gc.collect()

    def __iadd__(self, other):
        """Add a second matrix to this matrix."""
        # note: we don't check that the cells and features match,
        #       the other matrix must have a subset of self for this to be safe.
        #       This method is intended only to be used in the aggregate_hdfs
        #       method which first computes the full set of features and cells.
        feature_rows = self._s_features[
            np.searchsorted(self.features, other.features, sorter=self._s_features)]
        cell_cols = self._s_cells[
            np.searchsorted(self.cells, other.cells, sorter=self._s_cells)]

        if self.sparse:
            # convert other.matrix to COO if not already
            sub = other.matrix.tocoo()
            # remap row and col indices to self index space
            global_row = feature_rows[sub.row]
            global_col = cell_cols[sub.col]
            update = scipy.sparse.coo_matrix(
                (sub.data, (global_row, global_col)),
                shape=self.matrix.shape)

            self._matrix = self.matrix + update
        else:
            self.matrix[
                np.repeat(feature_rows, len(cell_cols)),
                np.tile(cell_cols, len(feature_rows))] += other.matrix.flatten()

        self._saturation += other.saturation
        return self

    def _remove_elements(self, feat_mask=None, cell_mask=None):
        """Remove matrix elements using masks and inplace copies."""
        if self.sparse:
            return self._remove_elements_sparse(feat_mask, cell_mask)

        # shuffle the rows (columns) we want to the front, then take a view.
        # Don't bother doing anything if the mask is everything
        if feat_mask is not None:
            sum_mask = np.count_nonzero(feat_mask)
            if sum_mask == 0:
                raise ValueError(
                    "All features would be removed, try altering filter thresholds.")
            if sum_mask != len(self.features):
                for i, x in enumerate(np.argwhere(feat_mask)):
                    self._matrix[i, :] = self._matrix[x, :]
                self._matrix = self._matrix[:i+1]
                self._features = self._features[feat_mask]
                self._s_features = np.argsort(self._features)

        if cell_mask is not None:
            sum_mask = np.count_nonzero(cell_mask)
            if sum_mask == 0:
                raise ValueError(
                    "All cells would be removed, try altering filter thresholds.")
            if sum_mask != len(self.cells):
                for j, x in enumerate(np.argwhere(cell_mask)):
                    self._matrix[:, j] = self.matrix[:, x].squeeze()
                self._matrix = self._matrix[:, :j+1]
                self._cells = self._cells[cell_mask]
                self._s_cells = np.argsort(self._cells)

        return self._matrix

    def _remove_elements_sparse(self, feat_mask=None, cell_mask=None):
        """Sparse-safe version of _remove_elements."""
        if feat_mask is not None:
            feat_mask = np.asarray(feat_mask).ravel()
            if feat_mask.ndim != 1 or feat_mask.shape[0] != len(self.features):
                raise ValueError("Invalid feature mask shape.")
            if not feat_mask.any():
                raise ValueError(
                    "All features would be removed, try altering filter thresholds.")

        if cell_mask is not None:
            cell_mask = np.asarray(cell_mask).ravel()
            if cell_mask.ndim != 1 or cell_mask.shape[0] != len(self.cells):
                raise ValueError("Invalid cell mask shape.")
            if not cell_mask.any():
                raise ValueError(
                    "All cells would be removed, try altering filter thresholds.")

        if feat_mask is not None and \
                np.count_nonzero(feat_mask) != len(self.features):
            self._matrix = self._matrix[feat_mask, :]
            self._features = self._features[feat_mask]
            self._s_features = np.argsort(self._features)

        if cell_mask is not None and \
                np.count_nonzero(cell_mask) != len(self.cells):
            self._matrix = self._matrix[:, cell_mask]
            self._cells = self._cells[cell_mask]
            self._s_cells = np.argsort(self._cells)

        return self._matrix

    def bin_cells_by_coordinates(self, bin_size=4, inplace=False):
        """Aggregate expression data across spatial bins."""
        if not self.is_visium_hd:
            raise ValueError(
                "This method is only applicable to Visium HD matrices with "
                "coordinates in the cell names.")
        self.logger.info("Binning cells by coordinates with bin size %d.", bin_size)
        self.logger.info(f"Cell barcodes look like: '{self._cells[0]}'.")

        pattern = re.compile(rb".*_(\d+)_(\d+)-1")
        match = pattern.match(self._cells[0])
        if not match:
            raise ValueError(f"Unexpected cell format: {self._cells[0]}")

        # extract coordinates, assume format at this point as we already
        # have checked one with regex
        coords = np.array([
            (int(cell[8:13]), int(cell[14:19]))
            for cell in self._cells
        ])
        binned_coords = coords // bin_size

        # group cell indices by (bx, by)
        bin_to_cells = defaultdict(list)
        for cell_idx, (bx, by) in enumerate(binned_coords):
            bin_to_cells[(bx, by)].append(cell_idx)

        # sort bin keys and assign column indices, this means the output
        # 'cell' list will always be a row major (in space)
        sorted_bins = sorted(bin_to_cells.keys())
        new_cell_names = [f"bin_{bx}_{by}".encode() for (bx, by) in sorted_bins]
        bin_key_to_column = {bin_key: i for i, bin_key in enumerate(sorted_bins)}

        # construct new matrix as COO. This is the fastest way to
        # aggregate the data. Where's the sum? Turns out that one of the
        # features of COO is that duplicate (row, col) pairs are summed
        # lazily when the matrix is operated on
        self.logger.info("Aggregating expression data across bins.")

        if self.sparse:
            # TODO: this may not be necessary, it was done because a previous
            #       implementation was slicing a lot along columns (e.g. selecting
            #       ranges of cells)
            self.logger.info("Converting to CSC")
            self._matrix = self.matrix.tocsc()
            self.logger.info("Matrix converted to CSC format.")

            row, col, data = [], [], []
            for bin_key, cell_indices in bin_to_cells.items():
                col_idx = bin_key_to_column[bin_key]
                for cell_idx in cell_indices:
                    start, end = (
                        self.matrix.indptr[cell_idx], self.matrix.indptr[cell_idx + 1])
                    rows = self.matrix.indices[start:end]
                    vals = self.matrix.data[start:end]
                    row.extend(rows)
                    col.extend([col_idx] * len(rows))
                    data.extend(vals)

            # pull together information into a sparse matrix
            binned_matrix = scipy.sparse.coo_matrix(
                (data, (row, col)),
                shape=(self.matrix.shape[0], len(sorted_bins))
            ).tocsr()
        else:
            # dense matrix, we can just sum the values
            binned_matrix = np.zeros(
                (self.matrix.shape[0], len(sorted_bins)), dtype=self._dtype)
            for bin_key, cell_indices in bin_to_cells.items():
                col_idx = bin_key_to_column[bin_key]
                binned_matrix[:, col_idx] = self.matrix[
                    :, cell_indices].sum(axis=1, dtype=self._dtype)

        self.logger.info("Finished binning cells by coordinates.")
        self.logger.info(
            f"Shape converted from {self._matrix.shape} to {binned_matrix.shape}")

        # just check these, because sparse matrices drive me crazy
        if binned_matrix.shape[0] != len(self.features):
            raise ValueError(
                f"Row mismatch: {binned_matrix.shape[0]} vs {len(self.features)}")
        if binned_matrix.shape[1] != len(new_cell_names):
            raise ValueError(
                f"Column mismatch: {binned_matrix.shape[1]} vs {len(new_cell_names)}")
        rtn = None
        if inplace:
            rtn = self
            self._matrix = binned_matrix
            self._cells = np.array(new_cell_names, dtype=bytes)
            # we've defined this to be true
            self._s_cells = np.arange(len(self._cells), dtype=int)
            # features remains the same
            self._fh = None
        else:
            rtn = ExpressionMatrix(
                matrix=binned_matrix,
                features=self.features,
                cells=np.array(new_cell_names, dtype=bytes),
                sparse=scipy.sparse.issparse(binned_matrix),
                dtype=self._dtype
            )
        return rtn


class Saturation:
    """Saturation calculator."""

    DTYPE = [("frac", float), ("n_reads", int), ("n_umis", int)]
    FRACTIONS = [
        1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1,
        0.05, 0.04, 0.03, 0.02, 0.01]

    def __init__(self, data):
        """Initialise saturation class."""
        if data.dtype != self.DTYPE:
            raise ValueError(
                f"Data must be a structured array with dtype {self.DTYPE}, "
                f"got {data.dtype}.")
        self._data = data

    @classmethod
    def from_tags(cls, df, random_state=42):
        """Sample per-read tags to derive no. UMIs vs no. reads."""
        arr = np.zeros(len(cls.FRACTIONS), dtype=cls.DTYPE)
        last_frac = None
        for i, frac in enumerate(cls.FRACTIONS):
            p = frac if last_frac is None else frac / last_frac
            df = df.sample(frac=p, random_state=random_state)
            arr["frac"][i] = frac
            arr["n_reads"][i] = len(df)
            # It's possible to have collisions of UMIs in a chunk that come from
            # the different genes and barcodes. This will be rare, but we should
            # still group to be correct.
            arr["n_umis"][i] = len(
                df.groupby(['corrected_umi', 'corrected_barcode', 'gene']))
            last_frac = frac
        return cls(arr)

    @classmethod
    def from_arr(cls, arr):
        """Initialize from structured array of pre-calculated data."""
        return cls(arr)

    @classmethod
    def zero(cls):
        """Create a zeroed saturation object."""
        sat = cls(np.zeros(len(cls.FRACTIONS), dtype=cls.DTYPE))
        sat._data['frac'] = cls.FRACTIONS
        return sat

    @property
    def saturation_results(self):
        """Get the saturation summary."""
        df = pd.DataFrame(self._data)
        df['saturation'] = 1 - (df['n_umis'] / df['n_reads'])
        df["saturation"].fillna(0, inplace=True)
        return df[['frac', 'n_umis', 'saturation', 'n_reads']]

    def __iadd__(self, other):
        """Add a second Saturation object."""
        if not np.allclose(self._data["frac"], other._data["frac"]):
            raise ValueError("Saturation fractions do not match.")
        self._data["n_reads"] += other._data["n_reads"]
        self._data["n_umis"] += other._data["n_umis"]
        return self
