#!/usr/bin/env python
"""Umap reduce."""
from pathlib import Path

import pandas as pd
from sklearn.decomposition import IncrementalPCA
import umap

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("umap_reduce")

    # Positional mandatory arguments
    parser.add_argument(
        "matrix",
        help="Gene expression matrix: rows=genes/transcripts, \
            columns=barcodes, values=UMIs",
    )
    parser.add_argument(
        "--pcn",
        help="Number of principal components to generate prior to UMAP \
        dimensionality reduction [100]",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--dimensions",
        help="Number of dimensions to reduce to [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--min_dist",
        help="Minimum distance apart that points are allowed to be in the 2D \
        map [0.1]",
        type=float,
        default=0.1,
    )

    parser.add_argument(
        "--n_neighbors",
        help="Number of neighbors to look at when learning the manifold \
        structure [15]",
        type=int,
        default=15,
    )

    parser.add_argument(
        "--output",
        help="UMAP TSV output file path.",
        type=Path
    )

    return parser


def run_pca(x, args):
    """Run PCA to generate args.pcn principal components."""
    # pcn must not exceed N samples or features
    pcn = min(args.pcn, x.shape[0], x.shape[1])
    transformer = IncrementalPCA(n_components=pcn)
    return transformer.fit_transform(x)


def run_umap(x_pca, index, args):
    """Create UMAP projections."""
    reducer = umap.UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        n_components=args.dimensions,
        verbose=2
    )

    # If there's only a single transcript column, the reducer
    # step will fail. For now just write an empty file.
    # TODO: Have a better way of detcting if transcript data does not have
    # enough transciript columns. Probably only an issue with test data
    try:
        x_embedded = reducer.fit_transform(x_pca)
    except TypeError:
        return None
    else:
        cols = [f"D{i+1}" for i in range(args.dimensions)]
        df_umap = pd.DataFrame(x_embedded, columns=cols, index=index)
        return df_umap


def main(args):
    """Run entry point."""
    logger = get_named_logger('UmapReduce')
    df = pd.read_csv(args.matrix, delimiter="\t", index_col=0)

    # Switch barcodes columns to rows
    x = df.transpose()

    logger.info(
        f"Running PCA: {x.shape[1]} features --> \
            {args.pcn} prinicipal components")
    x_pca = run_pca(x, args)

    logger.info(
        f"Running UMAP: {x_pca.shape[1]} features --> \
            {args.dimensions} dimensions")

    df_umap = run_umap(x_pca, x.index, args)
    if df_umap is not None:
        df_umap.to_csv(
            args.output, sep="\t", index=True, index_label="barcode")
    else:
        open(args.output, 'w').close()


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
