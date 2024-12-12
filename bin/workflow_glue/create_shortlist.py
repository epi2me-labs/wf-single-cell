"""Determine valid cell barcodes using frequency methods."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import functools
import itertools
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npm
import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101

logger = get_named_logger('CreatSList')


def argparser():
    """Create argument parser."""
    parser = wf_parser("create_shortlist")

    parser.add_argument(
        "input", type=Path,
        help="Directory of TSV tag files or counts store in pickle files.")

    parser.add_argument(
        "output", type=Path,
        help="Output filepath for barcode shortlist.")

    parser.add_argument(
        "summary", type=Path,
        help="Output filepath for TSV summary file.")

    grp = parser.add_mutually_exclusive_group(required=True)
    grp.add_argument(
        "--long_list", type=Path,
        help="Path to 10X whitelist ('long list').")
    grp.add_argument(
        "--counts", action='store_true',
        help="Input is pre-aggregated barcode counts stored as "
             "TSV files with columns 'barcode' and 'count'.")

    grp = parser.add_argument_group(title="Cell filtering")
    grp.add_argument(
        "--method", default="quantile",
        choices=["quantile", "fixed", "abundance", "distance"],
        help="Method for calculating knee position.")
    grp.add_argument(
        "--exp_cells", type=int, default=500,
        help="Used with 'quantile' method. An estimated cell count.")
    grp.add_argument(
        "--cell_count", type=int, default=500,
        help="Used with 'fixed' method. Select N most abundance cells.")
    grp.add_argument(
        "--read_count", type=int, default=10000,
        help="Used with 'abundance' method to pick cells using minimum read count.")

    grp = parser.add_argument_group("Optional outputs.")
    parser.add_argument(
        "--plot", type=Path,
        help="Knee plot filename")
    parser.add_argument(
        "--counts_out", type=Path,
        help="Barcode counts TSV file.")

    parser.add_argument(
        "--threads", type=int, default=4,
        help="Worker threasd for processing barcode data.")

    parser.add_argument(
        "--min_qv", type=int, default=15,
        help="Minimum base quality values of shortlisted exact hits to long_list.")

    parser.add_argument(
        "--no_cell_filter", action='store_true',
        help="Do not apply cell count thresholding to aggregated barcode counts.")

    return parser


def get_knee_quantile(count_array, exp_cells=500, percentile=0.95, fraction=20):
    """Quantile-based method for thresholding the cell barcode whitelist.

    This method is adapted from the following preprint:
    Yupei You, Yair D.J. Prawer, Ricardo De Paoli-Iseppi, Cameron P.J. Hunt,
    Clare L. Parish, Heejung Shim, Michael B. Clark. Identification of cell
    barcodes from long-read single-cell RNA-seq with BLAZE. biorxiv. 2022.
    doi: https://doi.org/10.1101/2022.08.16.504056
    """
    exp_cells = min(exp_cells, len(count_array))
    top_count = np.sort(count_array)[::-1][:exp_cells]
    threshold = np.quantile(top_count, percentile) / fraction
    return threshold


def get_knee_distance(values):
    """Get knee distance.

    This function is based on
    https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-
    point-on-a-curve and https://dataplatform.cloud.ibm.com/analytics
    /notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=
    6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1
    The idea is to draw a line from the first to last point on the
    cumulative counts curve and then find the point on the curve
    which is the maximum distance away from this line
    """
    # noqa
    # get coordinates of all the points
    n_points = len(values)
    all_coord = np.vstack((range(n_points), values)).T

    # get the first point
    first_point = all_coord[0]
    # get vector between first and last point - this is the line
    line_vec = all_coord[-1] - all_coord[0]
    line_vec_norm = line_vec / np.sqrt(np.sum(line_vec ** 2))

    # find the distance from each point to the line:
    # vector between all points and first point
    vec_from_first = all_coord - first_point

    # To calculate the distance to the line, we split vec_from_first into two
    # components, one that is parallel to the line and one that
    # is perpendicular.
    # Then, we take the norm of the part that is perpendicular to the line and
    # get the distance.
    # We find the vector parallel to the line by projecting vec_from_first onto
    # the line. The perpendicular vector is
    #   vec_from_first - vec_from_first_parallel
    # We project vec_from_first by taking the scalar product of the vector with
    # the unit vector that points in the direction of the line (this gives us
    # the length of the projection of vec_from_first onto the line). If we
    # multiply the scalar product by the unit vector, we have
    # vec_from_first_parallel

    scalar_product = np.sum(
        vec_from_first *
        npm.repmat(
            line_vec_norm,
            n_points,
            1),
        axis=1)
    vec_from_first_parallel = np.outer(scalar_product, line_vec_norm)
    vec_to_line = vec_from_first - vec_from_first_parallel

    # distance to line is the norm of vec_to_line
    dist_to_line = np.sqrt(np.sum(vec_to_line ** 2, axis=1))

    # knee/elbow is the point with max distance value
    idx = np.argmax(dist_to_line)
    return idx


def make_kneeplot(
        counts, n_cells, output, no_threshold_line=False, title=None, max_points=10000):
    """Make kneeplot.

    :param counts: sorted counts of cells.
    :param n_cells: selection point to demarcate.
    :param output: output filepath.
    :param max_points: maximum number of points to plot.
    :param no_threshold_line: If true do not plot vertical threshold (Visium data)
    """
    # plot is done on logscale from most to least abundant
    x = np.arange(0, len(counts))
    y = counts[::-1]

    x_label = "Barcode rank" if no_threshold_line else "Cell barcode rank"

    fig = plt.figure(figsize=[6, 6])
    ax1 = fig.add_subplot(111)
    ax1.scatter(x, y, color="k", alpha=0.1, s=5)
    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim([1, 100000])
    ax1.set_xlabel(x_label)
    ax1.set_ylabel("Read count")
    ymax = ax1.get_ylim()[1]
    if not no_threshold_line:
        ax1.vlines(n_cells, ymin=1, ymax=ymax, linestyle="--", color="k")
    ax1.set_title(f"Barcodes above cutoff: {n_cells}")
    fig.tight_layout()
    fig.savefig(output)
    return None


def _filter_barcodes(fnames, whitelist, min_qv=15):
    """Worker function for `filter_barcodes`."""
    counts = Counter()
    for fname in fnames:
        df = pd.read_csv(fname, sep='\t', usecols=['CR', 'CY'], dtype=str)
        # Keep only barcodes with 100% match in whitelist and all high quality bases
        df = df[df.CR.isin(whitelist)]
        df['min_q'] = df.CY.apply(lambda x: min(ord(y) - 33 for y in x))
        df = df[df.min_q >= min_qv]
        counts.update(df.CR)
    return counts


def filter_barcodes(whitelist, input_dir, min_qv=15, threads=4):
    """Count barcodes present in files."""
    fnames = list(Path(input_dir).iterdir())

    def _groups():
        it = [iter(fnames)] * (len(fnames) // threads)
        for item in itertools.zip_longest(*it):
            yield [x for x in item if x is not None]
    logger.info(f"Processing {len(fnames)} barcode files.")

    counts = None
    if len(fnames) == 1 or threads == 1:  # _groups() doesn't like a list of one
        counts = _filter_barcodes(fnames, whitelist, min_qv)
    else:
        counts = Counter()
        worker = functools.partial(_filter_barcodes, whitelist=whitelist, min_qv=min_qv)
        with ProcessPoolExecutor(max_workers=threads) as executor:
            for result in executor.map(worker, _groups()):
                counts.update(result)
    return counts


def aggregate_counts(input_dir):
    """Aggregate precomputed barcode counts either from pickle or TSV."""
    fnames = list(Path(input_dir).iterdir())

    counts = Counter()
    for fname in fnames:
        c = pd.read_csv(
            fname, delimiter="\t",
            dtype={"barcode": str, "count": int})
        counts.update(dict(zip(c["barcode"], c["count"])))
    return counts


def find_threshold(counts, method, exp_cells=500, cell_count=5000, read_count=10000):
    """Calculate a threshold count for selecting cells.

    :param: Series of sorted counts.
    """
    idx = len(counts)
    if method == "quantile":  # flames method
        logger.info(f"Using quantile method with {exp_cells} cells.")
        threshold = get_knee_quantile(counts, exp_cells)
        idx = np.searchsorted(counts, threshold)
    elif method == "distance":  # fancy method
        logger.info("Using extremal method.")
        idx = get_knee_distance(counts)
    elif method == "fixed":  # take a fixed number of cells
        logger.info(f"Using fixed method, taking {cell_count} cells.")
        idx = min(len(counts), cell_count)
    elif method == "abundance":  # threshold
        logger.info(f"Using abundance method with threshold {read_count} reads.")
        idx = np.searchsorted(counts, read_count)
    else:
        raise ValueError(f"Unknown cell selection method: '{method}'.")
    return idx


def main(args):
    """Run entry point."""
    bc_counts = None
    if args.counts:
        logger.info("Loading precomputed counts.")
        bc_counts = aggregate_counts(args.input)
    else:
        logger.info("Reading whitelist.")
        whitelist = pd.read_csv(
            args.long_list, header=None,
            names=["barcode"], dtype=str)["barcode"].values

        logger.info("Collecting high quality matches.")
        bc_counts = filter_barcodes(
            whitelist, args.input, min_qv=args.min_qv, threads=args.threads)
    if len(bc_counts) == 0:
        raise ValueError("No good seed barcodes found.")

    logger.info("Writing counts of perfect hits.")
    hq_bcs = pd.DataFrame(bc_counts.most_common(), columns=["barcode", "count"])
    if args.counts_out is not None:
        hq_bcs.to_csv(args.counts_out, sep='\t', index=False)

    logger.info(f"Finding count threshold from {len(hq_bcs)} good seeds.")
    hq_bcs = hq_bcs.iloc[::-1]  # we want them ascending for searching

    total_hq_reads = sum(hq_bcs["count"])
    if args.no_cell_filter:
        # Use all the high quality barcodes as our shortlist.
        # This is the case for Visium data
        # No thresholding applied
        threshold = 0
        threshold_index = 0
        hq_reads_remaining = total_hq_reads
        shortlist = hq_bcs["barcode"]
        n_cells = len(shortlist)
        logger.info(f"Outputting {n_cells} with no thresholding applied to reads.")
    else:
        threshold_index = find_threshold(
            hq_bcs["count"], args.method,
            exp_cells=args.exp_cells, cell_count=args.exp_cells,
            read_count=args.read_count)
        n_cells = len(hq_bcs) - threshold_index
        # The threshold (number of cells)
        threshold = hq_bcs["count"].iat[threshold_index]
        shortlist = hq_bcs["barcode"].iloc[threshold_index:]
        # High quality reads remaining after cell filtering
        hq_reads_remaining = sum(hq_bcs["count"].iloc[threshold_index:])
        logger.info(f"Outputting {n_cells} with more than {threshold} reads.")

    with open(args.output, "wt") as f:
        f.write("\n".join(shortlist[::-1]))  # most common first
        f.write("\n")

    if args.plot is not None:
        logger.info(f"Generating knee plot: {args.plot}")
        make_kneeplot(
            hq_bcs["count"], n_cells, args.plot,
            no_threshold_line=args.no_cell_filter)

    (
        pd.DataFrame.from_dict(
            dict(
                high_quality_reads=total_hq_reads,
                high_quality_barcodes=[len(hq_bcs)],
                n_cells_after_filtering=[n_cells],
                high_q_reads_in_cells=[hq_reads_remaining],
                fraction_of_high_q_reads_in_cells=[
                    hq_reads_remaining / total_hq_reads] if total_hq_reads > 0 else [0],
                cell_count_read_threshold=[threshold]))
        .transpose()
        .to_csv(args.summary, sep='\t', header=False)
    )
