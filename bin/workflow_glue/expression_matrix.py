"""Expression counts matrix construction."""
import gzip
import os

import h5py
import numpy as np
import pandas as pd
import polars as pl
import scipy.io
import scipy.sparse


class ExpressionMatrix:
    """Representation of expression matrices."""

    def __init__(
            self, matrix=None, features=None, cells=None,
            fname=None, cache=False, dtype=int):
        """Create a matrix.

        Do not use the constructor directly.
        """
        self._matrix = matrix
        self._features = features
        self._cells = cells

        if fname is not None:
            self._fh = h5py.File(fname, 'r')
        self._cache = cache
        self._dtype = dtype

        if features is not None and cells is not None and matrix is None:
            self._matrix = np.zeros((len(features), len(cells)), dtype=self._dtype)

        if features is not None:
            self._s_features = np.argsort(features)
        if cells is not None:
            self._s_cells = np.argsort(cells)

    def __exit__(self):
        """Cleanup."""
        if self._fh is not None:
            self._fh.close()

    @classmethod
    def from_hdf(cls, name, cache=True, dtype=int):
        """Load a matrix from HDF file."""
        ma = cls(fname=name, cache=cache, dtype=dtype)
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
        df.drop(df.index[df[feature] == '-'], inplace=True)
        df = (
            df.groupby([feature, 'corrected_barcode'], observed=True)
            .nunique()['corrected_umi']
            .reset_index()  # Get feature back to a column
            .pivot(index=feature, columns='corrected_barcode', values='corrected_umi')
        )
        df.fillna(0, inplace=True)
        return cls(
            df.to_numpy(dtype=int),
            df.index.to_numpy(dtype=bytes),
            df.columns.to_numpy(dtype=bytes),
            dtype=int)

    @classmethod
    def aggregate_hdfs(cls, fnames, dtype=int):
        """Aggregate a set of matrices stored in HDF."""
        if len(fnames) == 1:
            return cls.from_hdf(fnames[0], dtype=dtype)
        features = set()
        cells = set()
        for fname in fnames:
            ma = cls.from_hdf(fname)
            features.update(ma.features)
            cells.update(ma.cells)

        # sort by names, just to be nice, doesn't guarantee matrices can be compared
        full_matrix = cls(
            features=np.array(sorted(features), dtype=bytes),
            cells=np.array(sorted(cells), dtype=bytes),
            dtype=dtype)
        for fname in fnames:
            ma = cls.from_hdf(fname, dtype=dtype)
            full_matrix + ma
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
        elif self._fh is not None:
            features = self._fh['features'][()].astype(bytes)
            if self._cache:
                self._features = features
            return features
        else:
            raise RuntimeError("Matrix not initialized.")

    @property
    def cells(self):
        """An array of cell names."""
        if self._cells is not None:
            return self._cells
        elif self._fh is not None:
            cells = self._fh['cells'][()].astype(bytes)
            if self._cache:
                self._cells = cells
            return cells
        else:
            raise RuntimeError("Matrix not initialized.")

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
    def matrix(self):
        """Expression matrix array."""
        if self._matrix is not None:
            return self._matrix
        elif self._fh is not None:
            matrix = self._fh['matrix'].astype(self._dtype)[()]
            if self._cache:
                self._matrix = matrix
            return matrix
        else:
            raise RuntimeError("Matrix not initialized.")

    @property
    def mean_expression(self):
        """The mean expression of features per cell."""
        return self.matrix.mean(axis=0)

    def normalize(self, norm_count):
        """Normalize total cell weight to fixed cell count."""
        # cell_count / cell_total = X / <norm_count>
        total_per_cell = np.sum(self.matrix, axis=0, dtype=float)
        scaling = norm_count / total_per_cell
        if np.issubdtype(self._dtype, float):
            self._matrix *= scaling
        else:
            self._matrix = np.multiply(self.matrix, scaling, dtype=float)
        return self

    def remove_cells(self, threshold):
        """Remove cells with few features present."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_cells`.")
        n_features = np.count_nonzero(self.matrix, axis=0)
        mask = n_features >= threshold
        self._remove_elements(cell_mask=mask)
        return self

    def remove_features(self, threshold):
        """Remove features that are present in few cells."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_features`.")
        n_cells = np.count_nonzero(self.matrix, axis=1)
        mask = n_cells >= threshold
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
        n_features = np.count_nonzero(self.matrix, axis=0)
        cell_mask = n_features >= cell_thresh

        n_cells = np.count_nonzero(self.matrix, axis=1)
        feat_mask = n_cells >= feature_thresh

        self._remove_elements(feat_mask=feat_mask, cell_mask=cell_mask)

        return self

    def remove_skewed_cells(self, threshold, prefixes, fname=None, label=None):
        """Remove cells with overabundance of feature class."""
        if self.matrix.size == 0:
            raise ValueError("Matrix is zero-sized on entry to `remove_skewed_cells`.")
        sel = self.find_features(prefixes)
        sel_total = self.matrix[sel].sum(axis=0, dtype=float)
        total = self.matrix.sum(axis=0, dtype=float)
        sel_total /= total
        mask = sel_total <= threshold
        if fname is not None:
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
        np.log1p(self.matrix, out=self.matrix)
        self._matrix /= np.log(10)
        return self

    @staticmethod
    def write_matrix(fname, matrix, index, col_names, index_name="feature"):
        """Write a matrix (or vector) to TSV.

        :param fname: filename.
        :param matrix: matrix to write.
        :param index: index column.
        :param col_names: names for matrix columns.
        :param index_name: name for index column.
        """
        # ok, this is kinda horrible, but its fast
        # I don't think this is actually zero-copy because data is
        # C-contiguous not F-contiguous.
        df = pd.DataFrame(matrix, copy=False, columns=col_names, index=index)
        df.index.rename(index_name, inplace=True)
        df = pl.from_pandas(df, include_index=True)
        df.write_csv(fname, separator="\t")

    def __add__(self, other):
        """Add a second matrix to this matrix."""
        cols = self._s_features[
            np.searchsorted(self.features, other.features, sorter=self._s_features)]
        rows = self._s_cells[
            np.searchsorted(self.cells, other.cells, sorter=self._s_cells)]
        self.matrix[
            np.repeat(cols, len(rows)),
            np.tile(rows, len(cols))] += other.matrix.flatten()

    def _remove_elements(self, feat_mask=None, cell_mask=None):
        """Remove matrix elements using masks and inplace copies."""
        # shuffle the rows (columns) we want to the front, then take a view.
        # Don't bother doing anything if the mask is everything
        if feat_mask is not None:
            sum_mask = sum(feat_mask)
            if sum_mask == 0:
                raise ValueError(
                    "All features would be removed, try altering filter thresholds.")
            if sum_mask != len(self.features):
                for i, x in enumerate(np.argwhere(feat_mask)):
                    self._matrix[i, :] = self._matrix[x, :]
                self._matrix = self._matrix[:i+1]
                self._features = self._features[feat_mask]

        if cell_mask is not None:
            sum_mask = sum(cell_mask)
            if sum_mask == 0:
                raise ValueError(
                    "All cells would be removed, try altering filter thresholds.")
            if sum_mask != len(self.cells):
                for j, x in enumerate(np.argwhere(cell_mask)):
                    self._matrix[:, j] = self.matrix[:, x].squeeze()
                self._matrix = self._matrix[:, :j+1]
                self._cells = self._cells[cell_mask]

        return self._matrix
