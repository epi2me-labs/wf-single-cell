"""Umap reduce."""
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
import umap

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("umap_reduce")

    parser.add_argument(
        "matrix",
        help="Gene expression matrix: rows=genes/transcripts, "
             "columns=barcodes, values=UMIs")
    parser.add_argument(
        "output", type=Path,
        help="UMAP TSV output file path.")

    parser.add_argument(
        "--pcn", type=int, default=100,
        help="Number of principal components to generate prior to UMAP")

    parser.add_argument(
        "--dimensions", type=int, default=2,
        help="Number of dimensions in UMAP embedding")

    parser.add_argument(
        "--min_dist", type=float, default=0.1,
        help="Minimum distance parameter of UMAP")

    parser.add_argument(
        "--n_neighbors", type=int, default=15,
        help="Number of neighbors parameter of UMAP")

    return parser


def main(args):
    """Run entry point."""
    logger = get_named_logger('UmapReduce')

    # find the numpy of columns, since the numpy API doesn't allow skipping cols
    names = np.loadtxt(
        args.matrix, delimiter="\t", dtype=str, max_rows=1)
    names = names[1:]  # first is transcript/gene
    n_barcodes = len(names)

    logger.info("Expression matrix has {n_barcodes} cells.")
    logger.info("Reading entire matrix.")
    mat = np.loadtxt(
        args.matrix, delimiter="\t", dtype=float,
        skiprows=1, usecols=list(range(1, n_barcodes + 1)))
    mat = np.atleast_2d(mat).transpose()
    logger.info("Finished reading matrix.")

    logger.info(f"Expression matrix has shape: {mat.shape}")
    pcn = min(args.pcn, *mat.shape)
    model = PCA(n_components=pcn, copy=False)
    mat = model.fit_transform(mat)
    logger.info(f"PCA output matrix has shape: {mat.shape}")

    mapper = umap.UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        n_components=args.dimensions,
        verbose=0)
    embedding = mapper.fit_transform(mat)
    logger.info(f"UMAP Embedding has shape: {embedding.shape}")

    # would be nice to avoid a copy here, but the array is fairly small
    cols = [f"D{i+1}" for i in range(args.dimensions)]
    out = pd.DataFrame(embedding, columns=cols, index=names)
    out.to_csv(
        args.output, sep="\t", index=True, index_label="barcode")
