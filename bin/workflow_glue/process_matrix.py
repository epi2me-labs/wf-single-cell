"""Expression counts matrix construction."""
import argparse
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import umap

from .expression_matrix import ExpressionMatrix  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("exp_mat")

    parser.add_argument(
        "input", type=Path, nargs='+',
        help="TSV with read tag data or batched expression matrices in HDF.")
    parser.add_argument(
        "--feature", default="gene", choices=["gene", "transcript"],
        help="Feature to compute matrix. Only used when read tag input is given.")
    parser.add_argument(
        "--raw", default="raw_feature_bc_matrix",
        help="Output folder for raw counts MEX data.")
    parser.add_argument(
        "--processed", default="processed_feature_bc_matrix",
        help="Output folder for processed counts MEX data.")
    parser.add_argument(
        "--per_cell_expr", default="expression.mean-per-cell.tsv", type=Path,
        help="Output TSV for per-cell mean expression level.")
    parser.add_argument(
        "--per_cell_mito", default="expression.mito-per-cell.tsv", type=Path,
        help="Output TSV for per-cell mean expression level.")
    parser.add_argument(
        "--text", action="store_true", help=argparse.SUPPRESS)

    grp = parser.add_argument_group("Filtering")
    grp.add_argument(
        "--enable_filtering", action="store_true",
        help="Enable filtering of matrix.")
    grp.add_argument(
        "--min_features", type=int, default=100,
        help="Filter out cells that contain fewer features than this.")
    grp.add_argument(
        "--min_cells", type=int, default=3,
        help="Filter out features that are observed in fewer than this "
             "number of cells")
    grp.add_argument(
        "--max_mito", type=int, default=5,
        help="Filter out cells where more than this percentage of counts "
             "belong to mitochondrial features.")
    grp.add_argument(
        "--mito_prefixes", default=["MT-"], nargs='*',
        help="prefixes to identify mitochondrial features.")
    grp.add_argument(
        "--norm_count", type=int, default=10000,
        help="Normalize to this number of counts per cell as "
             "is performed in CellRanger.")
    grp.add_argument(
        "--filtered_mex", default="filtered_feature_bc_matrix",
        help="Output folder for raw counts MEX data.")

    grp = parser.add_argument_group("UMAP creation")
    grp.add_argument(
        "--enable_umap", action="store_true",
        help="Perform UMAP on matrix.")
    grp.add_argument(
        "--umap_tsv", default="expression.umap.tsv", type=Path,
        help=(
            "UMAP TSV output file path. If --replicates is greater than 1 "
            "files will be named: name.index.tsv."))
    grp.add_argument(
        "--replicates", type=int, default=1,
        help="Number of UMAP replicated to perform.")
    grp.add_argument(
        "--pcn", type=int, default=100,
        help="Number of principal components to generate prior to UMAP")
    grp.add_argument(
        "--dimensions", type=int, default=2,
        help="Number of dimensions in UMAP embedding")
    grp.add_argument(
        "--min_dist", type=float, default=0.1,
        help="Minimum distance parameter of UMAP")
    grp.add_argument(
        "--n_neighbors", type=int, default=15,
        help="Number of neighbors parameter of UMAP")

    return parser


def main(args):
    """Make feature x cell, UMI-deduplicated, counts matrix."""
    logger = get_named_logger('AggreMatrix')
    logger.info('Constructing count matrices')

    # converting to float on fly means we can save a copy when normalizing
    try:
        matrix = ExpressionMatrix.aggregate_tags(args.input, args.feature, dtype=float)
    except UnicodeDecodeError:
        matrix = ExpressionMatrix.aggregate_hdfs(args.input, dtype=float)

    logger.info("Removing unknown features.")
    if len(matrix.cells) == 0:
        raise ValueError("""The expression matrix contains no cells.
            This may indicate an issue with data quality or volume.
            Incorrectly specified 10x kits/versions and reference data can also lead to
            to removal of all data at this point.""")
    matrix.remove_unknown()

    logger.info("Writing raw counts to file.")
    if args.text:
        matrix.to_tsv(args.raw, args.feature)
    else:
        matrix.to_mex(args.raw, dtype=int)

    if args.enable_filtering:
        logger.info("Filtering, normalizing and log-transforming matrix.")
        matrix = (
            matrix
            .remove_cells_and_features(args.min_features, args.min_cells)
            .remove_skewed_cells(
                args.max_mito / 100, args.mito_prefixes,
                fname=args.per_cell_mito, label="mito_pct")
            .normalize(args.norm_count)
            .log_transform()
        )
        logger.info("Writing filtered matrix.")
        if args.text:
            matrix.to_tsv(args.processed, args.feature)
        else:
            matrix.to_mex(args.processed)
    else:
        logger.info("Normalizing and log-transforming matrix.")
        matrix.normalize(args.norm_count).log_transform()

    logger.info("Writing mean expression levels.")
    ExpressionMatrix.write_matrix(
        args.per_cell_expr,
        matrix.mean_expression, matrix.tcells, ['mean_expression'], index_name='CB')

    if args.enable_umap:
        # note, we're going to do things in place so ExpressionMatrix will
        # become modified (trimmed on feature axis, and transposed)
        mat = matrix._matrix
        pcn = min(args.pcn, *mat.shape)
        matrix._features = np.array([f"pca_{i}" for i in range(pcn)])
        logger.info(f"Performing PCA on matrix of shape: {matrix.matrix.shape}")
        model = PCA(n_components=pcn, copy=False)
        mat = model.fit_transform(mat.transpose())  # warning!
        logger.info(f"PCA output matrix has shape: {mat.shape}")

        for replicate in range(args.replicates):
            logger.info(f"Performing UMAP replicate {replicate + 1}.")
            mapper = umap.UMAP(
                n_neighbors=args.n_neighbors,
                min_dist=args.min_dist,
                n_components=args.dimensions,
                verbose=0)
            embedding = mapper.fit_transform(mat)
            logger.info(f"UMAP Embedding has shape: {embedding.shape}")

            # would be nice to avoid a copy here, but the array is fairly small
            fname = str(args.umap_tsv).replace('REPEAT', str(replicate))
            logger.info(f"Writing UMAP embedding {fname}.")
            cols = [f"D{i+1}" for i in range(args.dimensions)]
            out = pd.DataFrame(embedding, columns=cols, index=matrix.tcells)
            out.to_csv(fname, sep="\t", index=True, index_label="CB")
        matrix._matrix = mat.transpose()  # undo the PCA transpose

    logger.info("Done.")
