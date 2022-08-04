import argparse
import logging
import os
import sys

import matplotlib
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap, ListedColormap

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "umap",
        help="File containing the 2D UMAP projection of cell barcodes.",
        type=str,
    )

    parser.add_argument(
        "full_matrix",
        help="File containing the full expression matrix that was used for \
        the UMAP projection.",
        type=str,
        default=None,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output file name for UMAP plots [umap.png]",
        type=str,
        default="umap.png",
    )

    parser.add_argument(
        "-g",
        "--gene",
        help="Gene to annotate in UMAP plots (e.g. --gene=CD19). If not \
        specified, cells will be annotated with total UMI counts [None]",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--mito_genes",
        help="Annotate the expression level of all mitochondrial genes. This \
        is useful for identifying unreliable cells with high mitochondrial \
        gene expression. Overrides the --gene flag [False]",
        action="store_true",
    )

    # parser.add_argument(
    #     "-t",
    #     "--target_cells",
    #     help="List of cells to highlight in UMAP plots [None]",
    #     type=str,
    #     default=None,
    # )

    parser.add_argument(
        "-s", "--size", help="Size of markers [15]", type=int, default=15
    )

    parser.add_argument(
        "-a", "--alpha", help="Transpancy of markers [0.7]", type=float, default=0.7
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    if args.gene == "None":
        args.gene = None

    # Override --gene flag if --mito_genes is specified
    if args.mito_genes:
        args.gene = None

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def remove_top_right_axes(ax):
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.tick_params(axis="both", direction="in")
    ax.get_xaxis().tick_bottom()
    ax.get_yaxis().tick_left()


def scatterplot(df, values, args):
    """ """
    fig = plt.figure(figsize=[8, 8])
    ax = fig.add_axes([0.08, 0.08, 0.85, 0.85])

    # Create custom colormap from gray to blue
    colors = ["lightgray", "blue"]
    cmap_custom = LinearSegmentedColormap.from_list("bluegray", colors)

    if (values.name == "total") | (values.name == "mitochondrial"):
        cmap = cm.jet
    else:
        # cmap = cm.inferno_r
        cmap = cmap_custom

    plot = ax.scatter(
        df["D1"],
        df["D2"],
        s=args.size,
        edgecolor="None",
        c=values,
        cmap=cmap,
        alpha=args.alpha,
    )

    remove_top_right_axes(ax)

    plt.xlim([df["D1"].min() - 1, df["D1"].max() + 1])
    plt.ylim([df["D2"].min() - 1, df["D2"].max() + 1])

    n_cells = df.shape[0]
    if values.name == "highlight":
        title = f"{n_cells} cells: highlighted cells from {args.target_cells}"
    else:
        # title = f"{n_cells} cells: {values.name}"
        title = f"{values.name}"
        cbar = plt.colorbar(plot)
        cbar.set_label("Normalized expression", rotation=270, labelpad=15)

    ax.set_title(title)
    ax.set_xlabel("UMAP-1")
    ax.set_ylabel("UMAP-2")

    plt.savefig(args.output, dpi=300)


def get_expression(args):
    """ """
    df_f = (
        pd.read_csv(args.full_matrix, delimiter="\t")
        .rename(columns={"gene": "cell"})
        .set_index("cell")
    )
    df_f = df_f.transpose()

    # Create annotation dataframe with requested features
    df_annot = pd.DataFrame()
    df_annot["total"] = np.exp(df_f).sum(axis=1) - 1

    if args.gene:
        # Make sure requested gene is in the matrix
        if args.gene not in df_f.columns:
            logging.info(f"WARNING: gene {args.gene} not found in expression matrix!")
            fig = plt.figure(figsize=[8, 8])
            fig.add_axes([0.08, 0.08, 0.85, 0.85])
            plt.savefig(args.output)
            sys.exit()
        df_annot[args.gene] = df_f.loc[:, args.gene]
    elif args.mito_genes:
        # Make sure at least one mitochondrial gene is in the matrix
        mito_genes = [g for g in df_f.columns if g.startswith("MT-")]
        if len(mito_genes) == 0:
            logging.info("WARNING: No mitochondrial genes found in expression matrix!")
            fig = plt.figure(figsize=[8, 8])
            fig.add_axes([0.08, 0.08, 0.85, 0.85])
            plt.savefig(args.output)
            sys.exit()
        df_annot["mitochondrial"] = df_f.loc[:, mito_genes].mean(axis=1)
    df_annot.index.names = ["barcode"]
    return df_annot


def main(args):
    init_logger(args)

    df = pd.read_csv(args.umap, delimiter="\t").set_index("barcode")

    df_annot = get_expression(args)

    df = df.loc[df_annot.index]

    if (not args.gene) & (not args.mito_genes):
        logger.info("Plotting UMAP with total UMI counts")
        scatterplot(df, df_annot.loc[:, "total"], args)
    elif args.gene:
        logger.info(f"Plotting UMAP with {args.gene} annotations")
        scatterplot(df, df_annot.loc[:, args.gene], args)
    elif args.mito_genes:
        logger.info("Plotting UMAP with mitochrondrial gene annotations")
        scatterplot(df, df_annot.loc[:, "mitochondrial"], args)

    # if args.target_cells:
    #     logger.info(f"Plotting UMAP with highlighted cells from {args.target_cells}")
    #
    #     df_highlight = pd.read_csv(args.target_cells, header=None, names=["barcode"])
    #     df_highlight["barcode"] = df_highlight["barcode"].str.split(
    #         "-", n=0, expand=True
    #     )
    #     df_highlight["highlight"] = "red"
    #     df_highlight = df_highlight.set_index("barcode")
    #     df_highlight = pd.merge(df, df_highlight, on="barcode", how="left").fillna(
    #         "lightgray"
    #     )
    #     values = df_highlight.loc[:, "highlight"]
    #     scatterplot(df, values, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
