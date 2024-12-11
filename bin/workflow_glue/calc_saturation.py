"""Calculate saturation."""

from matplotlib import pyplot as plt
import pandas as pd
import polars as pl

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("Calculate satutation")

    parser.add_argument(
        "--read_tags",
        help="TSV file with read_id, gene, barcode, and UMI"
    )

    parser.add_argument(
        "--output",
        help="Output plot file with saturation curves [output.png]",
        default="output.png",
    )

    return parser


def plot_saturation_curves(res, args):
    """Plot saturation curves.

    Output a single file with two subplots:
    1. Median number of unique genes per cell
    2. Median number of unique UMIs per cell
    3. Sequencing saturation
    """
    fig = plt.figure(figsize=[15, 5])

    ax1 = fig.add_subplot(1, 3, 1)
    ax1.plot(
        res["reads_pc"],
        res["genes_pc"],
        marker=None,
        color="orange",
        linewidth=2)
    ax1.set_xlim([0, ax1.get_xlim()[1]])
    ax1.set_ylim([0, ax1.get_ylim()[1]])
    ax1.set_xlabel("Median reads per cell")
    ax1.set_ylabel("Genes per cell")
    ax1.set_title("Gene saturation")

    plt.setp(ax1.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax1.grid()

    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(
        res["reads_pc"],
        res["umis_pc"],
        marker=None,
        color="purple",
        linewidth=2)
    ax2.set_xlim([0, ax2.get_xlim()[1]])
    ax2.set_ylim([0, ax2.get_ylim()[1]])
    ax2.set_xlabel("Median reads per cell")
    ax2.set_ylabel("UMIs per cell")
    ax2.set_title("UMI saturation")
    plt.setp(ax2.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax2.grid()

    ax3 = fig.add_subplot(1, 3, 3)
    ax3.plot(
        res["reads_pc"],
        res["umi_sat"],
        marker=None,
        color="blue",
        linewidth=2)
    ax3.axhline(y=res["umi_sat"].max(), color="k", linestyle="--", linewidth=2)
    ax3.set_xlim([0, ax3.get_xlim()[1]])
    ax3.set_ylim([0, 1])
    ax3.set_xlabel("Median reads per cell")
    ax3.set_ylabel("Sequencing saturation")
    umi_sat = res.at[1, 'umi_sat']  # saturation at zero subsampling
    ax3.set_title(f"Sequencing saturation: {umi_sat:.2f}")
    plt.setp(ax3.xaxis.get_majorticklabels(), rotation=45, ha="right")
    ax3.grid()

    fig.tight_layout()
    fig.savefig(args.output)


def downsample_dataframe(df, fraction):
    """Downsample dataframe of read tags and tabulate genes and UMIs per cell."""
    logger = get_named_logger('ClcSat')

    logger.info(f"Doing {fraction}")
    df_scaled = df.sample(fraction=fraction)
    n_reads = df_scaled.shape[0]

    # Get the unique number of reads, genes and UMIs per cell barcode
    gb_cell = df_scaled.group_by("barcode")
    gb_cell_median = gb_cell.n_unique().median()
    genes_per_cell = gb_cell_median['gene'][0]
    umis_per_cell = gb_cell_median['umi'][0]
    # Since polars 0.20.5 groupby.count() has been renamed groupby.len()
    reads_per_cell = gb_cell.count().median()['count'][0]

    n_deduped_reads = df_scaled.group_by(['gene', 'barcode', 'umi']).count().shape[0]
    if n_reads < 1:
        umi_saturation = 0
    else:
        umi_saturation = 1 - (n_deduped_reads / n_reads)

    record = (
        (
            fraction,
            n_reads,
            reads_per_cell,
            genes_per_cell,
            umis_per_cell,
            umi_saturation,
        )
    )
    logger.info(f"Done saturation calculation for fraction {fraction}")
    return record


def run_jobs(args):
    """Create job to send off to workers, and collate results."""
    logger = get_named_logger('ClcSat')

    df = pl.read_csv(
        args.read_tags,
        separator='\t',
        columns=['corrected_barcode', 'corrected_umi', 'gene'],
        new_columns=['barcode', 'umi', 'gene'],
        low_memory=True,
        dtypes={
            'corrected_barcode': pl.Categorical,
            'corrected_umi': pl.Categorical,
            'gene': str}
    )

    df.filter((df['barcode'] != '-') & (df['umi'] != '-'))

    logger.info("Downsampling reads for saturation curves")
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
    for frac in fractions:
        records.append(downsample_dataframe(df, frac))

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
        index="downsamp_frac"
    )
    return res


def main(args):
    """Entry point."""
    result = run_jobs(args)
    plot_saturation_curves(result, args)
