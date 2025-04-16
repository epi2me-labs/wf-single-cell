"""Calculate saturation."""

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
        help="Output TSV file with saturation curves."
    )

    parser.add_argument(
        "--sample",
        help="sample ID/alias"
    )

    return parser


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

    res = pl.from_records(
        data=records,
        schema=[
            "downsamp_frac",
            "downsamp_reads",
            "reads_pc",
            "genes_pc",
            "umis_pc",
            "umi_sat",
        ]
    )
    res = res.with_columns(
        pl.lit(args.sample).alias("sample"),
    )
    res.write_csv(args.output, separator="\t")


def main(args):
    """Entry point."""
    run_jobs(args)
