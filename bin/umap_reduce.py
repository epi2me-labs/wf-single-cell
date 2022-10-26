#!/usr/bin/env python
"""Umap reduce."""
import argparse
import logging
from pathlib import Path

import pandas as pd
import umap


logger = logging.getLogger(__name__)


def parse_args():
    """Create argument parser."""
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "matrix",
        help="Gene expression matrix: rows=genes, columns=barcodes, \
        values=UMIs",
    )

    # Optional arguments
    parser.add_argument(
        "-o",
        "--output_prefix",
        help="write to current directory with this prefix",
        type=str,
        default="")

    parser.add_argument(
        "-d",
        "--dimensions",
        help="Number of dimensions to reduce to [2]",
        type=int,
        default=2,
    )

    parser.add_argument(
        "-f",
        "--feature_type",
        help="matrix is gene or transcript [gene]",
        default='gene',
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
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    parser.add_argument(
        "--num_umaps",
        help="Make multiple umap plots with different initial random states",
        type=int,
        default=10,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    """Initiate logger."""
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def main(args):
    """Run entry point."""
    init_logger(args)

    df = pd.read_csv(args.matrix, delimiter="\t").set_index(args.feature_type)

    # Switch from barcodes as columns to barcodes as rows
    X = df.transpose()

    logger.info(
        f"Running UMAP: {X.shape[1]} features --> \
            {args.dimensions} dimensions")

    for n in range(args.num_umaps):
        reducer = umap.UMAP(
            n_neighbors=args.n_neighbors,
            min_dist=args.min_dist,
            n_components=args.dimensions,
            verbose=2
        )
        outpath = Path() / f"{args.output_prefix}_{n}_umap.tsv"

        # For testing: If there's only a single transcript column the reducer
        # step will fail. For now just write an empty file.
        # TODO: Have a better way of detcting if transcript data does not have
        # enough transciript columns. Probably only an issue with test data
        try:
            X_embedded = reducer.fit_transform(X)
        except TypeError:
            open(outpath, 'w').close()
        else:
            cols = [f"D{i+1}" for i in range(args.dimensions)]
            df_umap = pd.DataFrame(X_embedded, columns=cols, index=X.index)

            df_umap.to_csv(
                outpath, sep="\t", index=True, index_label="barcode")


if __name__ == "__main__":
    args = parse_args()

    main(args)
