import argparse
import logging

import pandas as pd
from matplotlib import pyplot as plt
from tqdm import tqdm

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "gene_cell_umi",
        help="TSV file with read_id, gene, barcode, and UMI",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--output",
        help="Output plot file with saturation curves [output.png]",
        type=str,
        default="output.png",
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def plot_saturation_curves(res, umi_sat, args):
    """
    Output a single file with two subplots:
    1. Median number of unique genes per cell
    2. Median number of unique UMIs per cell
    """
    fig = plt.figure(figsize=[15, 5])

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(res["reads_pc"], res["genes_pc"], marker=None, color="orange", linewidth=2)
    ax1.set_xlim([0, ax1.get_xlim()[1]])
    ax1.set_ylim([0, ax1.get_ylim()[1]])
    ax1.set_xlabel("Median reads per cell")
    ax1.set_ylabel("Genes per cell")
    ax1.set_title("Gene saturation")

    # ax1.xaxis.set_major_locator(MultipleLocator(2000))
    # ax1.xaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax1.xaxis.set_minor_locator(MultipleLocator(500))

    # ax1.yaxis.set_major_locator(MultipleLocator(2000))
    # ax1.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    # ax1.yaxis.set_minor_locator(MultipleLocator(500))
    #
    # ax1.set_xlim([-5,20000])
    # ax1.set_ylim([-5,20000])
    # custom_legend(ax1)
    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax1.grid()

    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(res["reads_pc"], res["umis_pc"], marker=None, color="purple", linewidth=2)
    ax2.set_xlim([0, ax2.get_xlim()[1]])
    ax2.set_ylim([0, ax2.get_ylim()[1]])
    ax2.set_xlabel("Median reads per cell")
    ax2.set_ylabel("UMIs per cell")
    ax2.set_title("UMI saturation")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax2.grid()

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(res["reads_pc"], res["umi_sat"], marker=None, color="blue", linewidth=2)
    ax3.axhline(y=res["umi_sat"].max(), color="k", linestyle="--", linewidth=2)
    ax3.set_xlim([0, ax3.get_xlim()[1]])
    ax3.set_ylim([0, 1])
    ax3.set_xlabel("Median reads per cell")
    ax3.set_ylabel("Sequencing saturation")
    ax3.set_title(f"Sequencing saturation: {umi_sat:.2f}")
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax3.grid()

    fig.tight_layout()
    fig.savefig(args.output)


def calc_umi_saturation(df):
    """
    Sequencing Saturation = 1 - (n_deduped_reads / n_reads)
    """
    n_reads = df.shape[0]
    df["gene_bc_umi"] = df["gene"] + "_" + df["barcode"] + "_" + df["umi"]
    n_deduped_reads = df["gene_bc_umi"].nunique()
    saturation = 1 - (n_deduped_reads / n_reads)

    return saturation


def downsample_reads(df):
    """
    Downsample dataframe of reads and tabulate genes and UMIs per cell
    """
    fractions = [
        0.01,
        0.02,
        0.03,
        0.04,
        0.05,
        0.1,
        0.2,
        0.3,
        0.4,
        0.5,
        0.6,
        0.7,
        0.8,
        0.9,
        1.0,
    ]

    records = [(0.0, 0, 0, 0, 0, 0.0)]
    for fraction in tqdm(fractions, total=len(fractions)):
        df_ = df.sample(frac=fraction)
        downsamp_reads = df_.shape[0]
        # Get the unique number of reads, genes and UMIs per cell barcode
        reads_per_cell = df_.groupby("barcode")["read_id"].nunique().median()
        genes_per_cell = df_.groupby("barcode")["gene"].nunique().median()
        umis_per_cell = df_.groupby("barcode")["umi"].nunique().median()
        umi_sat = calc_umi_saturation(df_)
        records.append(
            (
                fraction,
                downsamp_reads,
                reads_per_cell,
                genes_per_cell,
                umis_per_cell,
                umi_sat,
            )
        )

    res = pd.DataFrame.from_records(
        records,
        columns=[
            "downsamp_frac",
            "downsamp_reads",
            "reads_pc",
            "genes_pc",
            "umis_pc",
            "umi_sat",
        ],
    )
    return res


def main(args):
    init_logger(args)

    df = pd.read_csv(args.gene_cell_umi, sep="\t")

    umi_sat = calc_umi_saturation(df)
    logger.info(f"Sequencing saturation: {umi_sat:.2f}")

    logger.info("Downsampling reads for saturation curves")
    res = downsample_reads(df)

    plot_saturation_curves(res, umi_sat, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
