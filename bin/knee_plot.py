import argparse
import logging
import operator
import sys
from collections import Counter

import matplotlib.pyplot as plt
import numpy as np
import numpy.matlib as npm
from scipy.signal import argrelextrema
from scipy.stats import gaussian_kde

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "barcodes",
        help="TSV file containing counts for uncorrected barcodes that are \
            present in the barcode superlist.",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "--knee_method",
        help="Method (distance/density) to use for calculating \
                        knee position [distance]",
        type=str,
        default="distance",
    )

    parser.add_argument(
        "--cell_count",
        help="Instead of empirically setting the cell count \
                        threshold from the knee plot using distance or density \
                        algorithms, use the top <N> cells based on read count. \
                        This option overrides the <knee_method> for generating \
                        the whitelist [None]",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--read_count_threshold",
        help="Instead of empirically setting the cell count \
                        threshold from the knee plot using distance or density \
                        algorithms, pick cells based on an explicit minimum \
                        read count per cell. This option overrides the \
                        <knee_method> for generating the whitelist [None]",
        type=int,
        default=None,
    )

    parser.add_argument(
        "--output_plot",
        help="Knee plot filename [kneeplot.png]",
        type=str,
        default="kneeplot.png",
    )

    parser.add_argument(
        "--output_whitelist",
        help="Barcode whitelist filename [ont_barcodes.tsv]",
        type=str,
        default="ont_barcodes.tsv",
    )

    parser.add_argument(
        "--ilmn_barcodes",
        help="Illumina barcodes filename, to be overlayed onto \
                        the knee plot of ONT barcodes [None]",
        type=str,
        default=None,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()
    init_logger(args)

    if (args.cell_count is not None) and (args.read_count_threshold is not None):
        raise Exception(
            "Cannot specify BOTH --cell_count and --read_count_threshold. Please pick one method."
        )

    if args.cell_count is not None:
        logging.info(
            f"Using explicit cell count (N = {args.cell_count}) for creating whitelist, not the distance or density algorithm"
        )

    if args.read_count_threshold is not None:
        logging.info(
            f"Using explicit reads per cell threshold (>= {args.read_count_threshold}) for creating whitelist, not the distance or density algorithm"
        )

    return args


def getKneeDistance(values):
    """
    This function is based on
    https://stackoverflow.com/questions/2018178/finding-the-best-trade-off-point-on-a-curve
    and https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1
    The idea is to draw a line from the first to last point on the
    cumulative counts curve and then find the point on the curve
    which is the maximum distance away from this line
    """

    # get coordinates of all the points
    nPoints = len(values)
    allCoord = np.vstack((range(nPoints), values)).T

    # get the first point
    firstPoint = allCoord[0]
    # get vector between first and last point - this is the line
    lineVec = allCoord[-1] - allCoord[0]
    lineVecNorm = lineVec / np.sqrt(np.sum(lineVec ** 2))

    # find the distance from each point to the line:
    # vector between all points and first point
    vecFromFirst = allCoord - firstPoint

    # To calculate the distance to the line, we split vecFromFirst into two
    # components, one that is parallel to the line and one that is perpendicular
    # Then, we take the norm of the part that is perpendicular to the line and
    # get the distance.
    # We find the vector parallel to the line by projecting vecFromFirst onto
    # the line. The perpendicular vector is vecFromFirst - vecFromFirstParallel
    # We project vecFromFirst by taking the scalar product of the vector with
    # the unit vector that points in the direction of the line (this gives us
    # the length of the projection of vecFromFirst onto the line). If we
    # multiply the scalar product by the unit vector, we have vecFromFirstParallel

    scalarProduct = np.sum(vecFromFirst * npm.repmat(lineVecNorm, nPoints, 1), axis=1)
    vecFromFirstParallel = np.outer(scalarProduct, lineVecNorm)
    vecToLine = vecFromFirst - vecFromFirstParallel

    # distance to line is the norm of vecToLine
    distToLine = np.sqrt(np.sum(vecToLine ** 2, axis=1))

    # knee/elbow is the point with max distance value
    idxOfBestPoint = np.argmax(distToLine)

    return distToLine, idxOfBestPoint


def getKneeEstimateDensity(
    cell_barcode_counts, expect_cells=False, cell_number=False, plotfile_prefix=None
):
    """estimate the number of "true" cell barcodes using a gaussian
    density-based method
    input:
         cell_barcode_counts = dict(key = barcode, value = count)
         expect_cells (optional) = define the expected number of cells
         cell_number (optional) = define number of cell barcodes to accept
         plotfile_prefix = (optional) prefix for plots
    returns:
         List of true barcodes
    """

    # very low abundance cell barcodes are filtered out (< 0.001 *
    # the most abundant)
    threshold = 0.001 * cell_barcode_counts.most_common(1)[0][1]

    counts = sorted(cell_barcode_counts.values(), reverse=True)
    counts_thresh = [x for x in counts if x > threshold]
    log_counts = np.log10(counts_thresh)

    # guassian density with hardcoded bw
    density = gaussian_kde(log_counts, bw_method=0.1)

    xx_values = 10000  # how many x values for density plot
    xx = np.linspace(log_counts.min(), log_counts.max(), xx_values)

    local_min = None

    if cell_number:  # we have a prior hard expectation on the number of cells
        threshold = counts[cell_number]

    else:
        local_mins = argrelextrema(density(xx), np.less)[0]
        local_mins_counts = []

        for poss_local_min in local_mins[::-1]:

            passing_threshold = sum(
                [
                    y > np.power(10, xx[poss_local_min])
                    for x, y in cell_barcode_counts.items()
                ]
            )
            local_mins_counts.append(passing_threshold)

            if not local_min:  # if we have selected a local min yet
                if expect_cells:  # we have a "soft" expectation
                    if (
                        passing_threshold > expect_cells * 0.1
                        and passing_threshold <= expect_cells
                    ):
                        local_min = poss_local_min

                else:  # we have no prior expectation
                    # TS: In abscence of any expectation (either hard or soft),
                    # this set of heuristic thresholds are used to decide
                    # which local minimum to select.
                    # This is very unlikely to be the best way to achieve this!
                    if poss_local_min >= 0.2 * xx_values and (
                        log_counts.max() - xx[poss_local_min] > 0.5
                        or xx[poss_local_min] < log_counts.max() / 2
                    ):
                        local_min = poss_local_min

        if local_min is not None:
            threshold = np.power(10, xx[local_min])

    if cell_number or local_min is not None:
        final_barcodes = set(
            [x for x, y in cell_barcode_counts.items() if y > threshold]
        )
    else:
        final_barcodes = None

    return final_barcodes, threshold


def get_barcode_counts(barcodes):
    barcode_counts = {}
    for line in open(barcodes, "r"):
        line = line.strip()
        seq = line.split("\t")[0]
        n = int(line.split("\t")[1])
        barcode_counts[seq] = n
    return barcode_counts


def apply_bc_cutoff(ont_bc_sorted, idx):
    return set(list(ont_bc_sorted.keys())[: idx + 1])


def write_ont_barcodes(cutoff_ont_bcs, args):
    with open(args.output_whitelist, "wt") as f:
        f.write("\n".join(list(cutoff_ont_bcs)))
        f.write("\n")


def make_kneeplot(ont_bc, ilmn_bc, conserved_bc, args):
    ont_bc_sorted = dict(
        sorted(ont_bc.items(), key=operator.itemgetter(1), reverse=True)
    )

    if args.ilmn_barcodes is not None:
        fig = plt.figure(figsize=[6, 10])
        ax1 = fig.add_subplot(311)
    else:
        fig = plt.figure(figsize=[6, 6])
        ax1 = fig.add_subplot(111)

    # Only plot 50 cells for a given number of reads. This dramatically reduces
    # the number of points to plot in the long tail, which is helpful when
    # making vectorized images.
    X = []
    Y = []
    n_counter = Counter()
    for i, (b, n) in enumerate(ont_bc_sorted.items()):
        n_counter[n] += 1
        if n_counter[n] <= 50:
            X.append(i)
            Y.append(n)

    ax1.scatter(
        X,
        Y,
        color="k",
        label="ONT",
        alpha=0.1,
        s=5,
    )

    ax1.set_xscale("log")
    ax1.set_yscale("log")
    ax1.set_xlim([1, 100000])
    ax1.set_xlabel("Cells")
    ax1.set_ylabel("Read count")
    ymax = ax1.get_ylim()[1]

    # Calculate the knee index
    if (args.cell_count is None) and (args.read_count_threshold is None):
        if args.knee_method == "distance":
            ont_counts = list(ont_bc_sorted.values())
            distToLine, idxOfBestPoint = getKneeDistance(ont_counts)
            cutoff_ont_bcs = apply_bc_cutoff(ont_bc_sorted, idxOfBestPoint)
        elif args.knee_method == "density":
            cutoff_ont_bcs, threshold = getKneeEstimateDensity(Counter(ont_bc))
            idxOfBestPoint, countOfBestPoint = min(
                enumerate(ont_bc_sorted.values()), key=lambda x: abs(x[1] - threshold)
            )
        else:
            print("Invalid value for --knee_method (distance, density)")
            sys.exit()
        logger.info(
            f"Writing {len(cutoff_ont_bcs)} cells to {args.output_whitelist} based on the {args.knee_method} algorithm"
        )
    elif args.cell_count is not None:
        cell_count = min(len(ont_bc_sorted), args.cell_count)
        cutoff_ont_bcs = set(list(ont_bc_sorted.keys())[:cell_count])
        idxOfBestPoint = cell_count
        logger.info(f"Writing top {cell_count} cells to {args.output_whitelist}")
    elif args.read_count_threshold is not None:
        cutoff_ont_bcs = set(
            [bc for bc, n in ont_bc_sorted.items() if n >= args.read_count_threshold]
        )
        idxOfBestPoint = len(cutoff_ont_bcs)
        logger.info(
            f"Writing {len(cutoff_ont_bcs)} cells with >= {args.read_count_threshold} reads to {args.output_whitelist}"
        )

    write_ont_barcodes(cutoff_ont_bcs, args)

    ax1.vlines(idxOfBestPoint, ymin=1, ymax=ymax, linestyle="--", color="k")
    ax1.set_title("Found {} cells using ONT barcodes".format(idxOfBestPoint))

    if args.ilmn_barcodes is not None:
        pct_ilmn_in_ont = 100 * len(cutoff_ont_bcs & ilmn_bc) / len(ilmn_bc)

        ax2 = fig.add_subplot(312)
        ax2.scatter(
            range(len(ont_bc.keys())),
            ont_bc_sorted.values(),
            color="k",
            label="ONT",
            alpha=0.1,
            s=5,
        )
        conserved_to_plot = []
        for i, (bc, n) in enumerate(ont_bc_sorted.items()):
            if conserved_bc.get(bc) is not None:
                conserved_to_plot.append((i, conserved_bc.get(bc)))
        x, y = zip(*conserved_to_plot)
        ax2.scatter(
            x, y, color="r", alpha=0.6, marker="+", s=15, label="Both ONT + ILMN"
        )
        ax2.vlines(idxOfBestPoint, ymin=1, ymax=ymax, linestyle="--", color="k")
        ax2.set_xscale("log")
        ax2.set_yscale("log")
        ax2.set_xlim([1, 100000])
        ax2.set_ylim([1, ymax])
        ax2.legend()
        ax2.set_xlabel("Cells")
        ax2.set_ylabel("Read count")
        ax2.set_title(
            "{:.1f}% of {} ILMN barcodes in {} ONT barcodes".format(
                pct_ilmn_in_ont, len(ilmn_bc), len(cutoff_ont_bcs)
            )
        )

        ax3 = fig.add_subplot(313)
        ont_only_to_plot = []
        ont_only_bcs = set()
        for i, (bc, n) in enumerate(ont_bc_sorted.items()):
            if conserved_bc.get(bc) is None:
                ont_only_to_plot.append((i, ont_bc_sorted.get(bc)))
                if bc in cutoff_ont_bcs:
                    ont_only_bcs.add(bc)
        x, y = zip(*ont_only_to_plot)
        ax3.scatter(x, y, color="b", label="ONT only (not in ILMN)", alpha=0.1, s=5)
        ax3.vlines(idxOfBestPoint, ymin=1, ymax=ymax, linestyle="--", color="k")
        ax3.set_xscale("log")
        ax3.set_yscale("log")
        ax3.set_xlim([1, 100000])
        ax3.set_ylim([1, ymax])
        ax3.legend()
        ax3.set_xlabel("Cells")
        ax3.set_ylabel("Read count")
        ax3.set_title(
            "{:.1f}% of {} ONT barcodes not in {} ILMN barcodes".format(
                100 * len(ont_only_bcs) / len(cutoff_ont_bcs),
                len(cutoff_ont_bcs),
                len(ilmn_bc),
            )
        )

    ax1.legend()
    fig.tight_layout()
    fig.savefig(args.output_plot)


def read_ilmn_barcodes(barcodes):
    bc = set([x.split("-")[0] for x in np.loadtxt(barcodes, dtype="str")])
    return bc


def intersect_ont_ilmn(ont_bc, ilmn_bc):
    conserved_bc = {}
    for bc, n in ont_bc.items():
        if bc in ilmn_bc:
            conserved_bc[bc] = n
    return conserved_bc


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def main(args):
    ont_bc = get_barcode_counts(args.barcodes)
    if args.ilmn_barcodes is not None:
        ilmn_bc = read_ilmn_barcodes(args.ilmn_barcodes)
        conserved_bc = intersect_ont_ilmn(ont_bc, ilmn_bc)
    else:
        ilmn_bc = {}
        conserved_bc = {}

    logger.info(f"Generating knee plot: {args.output_plot}")
    make_kneeplot(ont_bc, ilmn_bc, conserved_bc, args)


if __name__ == "__main__":
    args = parse_args()

    main(args)
