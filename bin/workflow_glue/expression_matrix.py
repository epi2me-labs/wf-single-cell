#!/usr/bin/env python
"""Expression counts matrix construction."""
from pathlib import Path

import numpy as np
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("exp_mat")

    parser.add_argument(
        "read_tags",
        help="TSV with read_id and associated tags.",
        type=Path
    )

    parser.add_argument(
        "feature_type",
        help="Feature to process.",
        choices=('gene', 'transcript')
    )

    parser.add_argument(
        "output_prefix",
        help="Output file prefix"
    )

    parser.add_argument(
        "--min_features",
        help="Filter out cells that contain fewer features than this.",
        type=int,
        default=100,
    )

    parser.add_argument(
        "--min_cells",
        help="Filter out features that are observed in fewer than this "
             "number of cells",
        type=int,
        default=3,
    )

    parser.add_argument(
        "--max_mito",
        help="Filter out cells where more than this percentage of counts "
             "belong to mitochondrial features.",
        type=int,
        default=5,
    )

    parser.add_argument(
        "--mito_prefixes",
        help="prefixes to identify mitochondrial features.",
        default=["MT-"],
        nargs='*'
    )

    parser.add_argument(
        "--norm_count",
        help="Normalize to this number of counts per cell as "
             "is performed in CellRanger.",
        type=int,
        default=10000,
    )

    return parser


def filter_cells(
        df, feature_type, min_features, max_mito, mito_prefixes, output_prefix):
    """Remove cells with small number of unique features."""
    logger = get_named_logger('FeatExpr')
    df = df.transpose()

    # Remove cells that have fewer than <min_features> unique features
    n_features = np.count_nonzero(df, axis=1)
    n_cells_prefilter = len(df)
    df = df[n_features >= min_features]
    n_dropped = n_cells_prefilter - len(df)
    logger.info(f"Dropping {n_dropped} cells with < {min_features} {feature_type}")

    # Filter cells with large numbers of mitochondrial features.
    # Note this does not apply to isoforms
    if feature_type == 'gene':
        mito_features = []
        for feature in df.columns:
            for mpre in mito_prefixes:
                if feature.find(mpre) == 0:
                    mito_features.append(feature)

        mito_total = df.loc[:, mito_features].sum(axis=1)
        mito_pct = 100 * mito_total / df["counts_per_cell"]
        mito_pct.to_csv(
            f"{output_prefix}_expression.mito.tsv", sep="\t")
        n_mito = mito_pct[mito_pct > max_mito].shape[0]
        logger.info(
            f"Dropping {n_mito} cells with > {max_mito}% mitochondrial reads")
        df = df[mito_pct <= max_mito]

    df = df.transpose()

    if df.shape[1] == 0:
        raise Exception("All cells have been filtered out!")

    return df


def normalize(df, norm_count):
    """Normalize."""
    logger = get_named_logger('FeatExpr')
    df = df.transpose()

    # cell_count / cell_total = X / <norm_count>
    logger.info(f"Normalizing counts to {norm_count} reads per cell")
    features = df.columns != "counts_per_cell"
    scaling = norm_count / df["counts_per_cell"]
    df.loc[:, features] = df.loc[:, features].mul(scaling, axis=0)
    df = df.transpose()

    return df


def filter_features(df, min_cells):
    """Filter features that are present in few cells."""
    logger = get_named_logger('FeatExpr')
    n_cells = np.count_nonzero(df, axis=1)
    n_features_prefilter = len(df)
    df = df[n_cells >= min_cells]

    if df.shape[0] == 0:
        raise Exception("All features have been filtered out!")

    n_features_dropped = n_features_prefilter - len(df)
    logger.info(
        f"Dropping {n_features_dropped} features observed in < {min_cells} cells")

    return df


def main(args):
    """
    Make feature x cell, UMI-deduplicated, counts matrix.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    """
    logger = get_named_logger('FeatExpr')
    logger.info('Constructing count matrices')

    df = pd.read_csv(
        args.read_tags,
        index_col=None,
        usecols=['CB', 'UB', args.feature_type],
        sep='\t',
        dtype='category'
    )

    # remove reads with no feature assigned
    df = df[df[args.feature_type] != '-']
    # Create the feature x cell matrix
    df = (
        df.groupby([args.feature_type, 'CB'], observed=True).nunique()['UB']
        .reset_index()  # Get feature back to a column
        .pivot(index=args.feature_type, columns='CB', values='UB')
    )
    df.fillna(0, inplace=True)

    # Write the raw counts matrix
    df.to_csv(
        f'{args.output_prefix}_expression.count.tsv',
        sep="\t"
    )

    # Process the raw counts matrix. Whether we filter cells or features first will have
    # an effect on the results. There's no specific reason for choosing cell filtering
    # first here.
    df.loc["counts_per_cell"] = df.sum(axis=0)
    df_proc = filter_cells(
        df, args.feature_type, args.min_features, args.max_mito,
        args.mito_prefixes, args.output_prefix)
    df_proc = filter_features(df_proc, args.min_cells)
    df_proc = normalize(df_proc, args.norm_count)
    df_proc.drop(["counts_per_cell"], axis=0, inplace=True)
    df_proc = np.log1p(df_proc) / np.log(10)

    # Write processed matrix
    df_proc.to_csv(
        f'{args.output_prefix}_expression.processed.tsv',
        sep="\t")


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
