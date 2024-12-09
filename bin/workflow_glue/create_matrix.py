"""Correct UMIs using clustering."""

import collections
import itertools
from pathlib import Path

from editdistance import eval as edit_distance
import numpy as np
import pandas as pd
from umi_tools import UMIClusterer

from .expression_matrix import ExpressionMatrix  # noqa: ABS101
from .tag_bam import BAM_TAGS  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101
from .sc_util import StatsSummary  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("cluster_umis")

    parser.add_argument(
        "chrom",
        help="Chromosome name")
    parser.add_argument(
        "barcode_tags", type=Path,
        help="Read tags TSV file.")
    parser.add_argument(
        "features", type=Path,
        help="TSV read gene/transcript assignments file.")
    grp = parser.add_argument_group("Output")
    grp.add_argument(
        "--tsv_out", type=Path,
        help="Output TSV containing a subset of read-tags in human-readable form")
    grp.add_argument(
        "--hdf_out", type=Path,
        help="Output filename for HDF matrix output. \
            Two files will be produced as filename.{gene, transcript}.ext")
    grp.add_argument(
        "--stats", type=Path,
        help="Output filename for JSON statistics summary. \
            Two files will be produced as filename.{gene, transcript}.ext")
    parser.add_argument(
        "--ref_interval", type=int, default=1000,
        help="Size of genomic window (bp) to assign as gene name if no gene \
            assigned by featureCounts.")

    return parser


def get_adj_list_directional_lev(self, umis, counts, threshold=2):
    """Use Levenshtein distance for UMIclustering instead of hamming.

    This function is to monkey-patch UMIClusterer._get_adj_list_directional
    """
    adj_list = {umi: [] for umi in umis}
    iter_umi_pairs = itertools.combinations(umis, 2)
    for umi1, umi2 in iter_umi_pairs:
        if edit_distance(umi1, umi2) <= threshold:
            if counts[umi1] >= (counts[umi2] * 2) - 1:
                adj_list[umi1].append(umi2)
            if counts[umi2] >= (counts[umi1] * 2) - 1:
                adj_list[umi2].append(umi1)

    return adj_list


def umi_clusterer_call(self, umis, threshold):
    """To replace UMIClusterer.__call__.

    Use this method to monkey-patch the UMICluterer.__call__ in order to remove the
    nessesity for all UMIs to be the same length, allowing for deletions in the UMIs.

    https://github.com/CGATOxford/UMI-tools/blob/c3ead0792ad590822ca72239ef01b8e559802d
    a9/umi_tools/network.py#L350
    """
    counts = umis
    umis = list(umis.keys())

    self.positions += 1
    number_of_umis = len(umis)
    self.total_umis_per_position += number_of_umis

    if number_of_umis > self.max_umis_per_position:
        self.max_umis_per_position = number_of_umis

    adj_list = self.get_adj_list(umis, counts, threshold)
    clusters = self.get_connected_components(umis, adj_list, counts)
    final_umis = [list(x) for x in self.get_groups(clusters, adj_list, counts)]

    return final_umis


# Monkey-patch the umi-tools clusterer with a modified method
# using Levenshtein instead of Hamming distance.
UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
UMIClusterer.__call__ = umi_clusterer_call


def cluster(umis):
    """Cluster UMIs.

    Search for UMI clusters within subsets of reads sharing the same corrected barcode
    and gene. In this way the search space is dramatically reduced.

    We are using the UMI-tools directional deduplication method (modified to use
    Levenshtein distance). Connections between nodes within a cluster are generated
    based on edit distance threshold and whether node A counts >= (2* node B counts).
    https://umi-tools.readthedocs.io/en/latest/the_methods.html

    """
    if len(umis) == 1:  # early return
        return umis
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)

    if len(clusters) == len(umis):  # no corrections
        return umis

    # create list of corrections
    umi_map = dict()
    for clust in clusters:
        if len(clust) > 1:
            for umi in clust[1:]:
                umi_map[umi] = clust[0]
    if len(umi_map) > 0:  # pd.Series.replace is weird slow
        umis = umis.replace(umi_map)
    return umis


def create_region_name(row, ref_interval):
    """Create a fake gene name from alignment coordinates."""
    # The idea here is to slice the reference into a grid and label reads
    # with the chunk that they overlap. Reads intersecting the same chunk
    # are then grouped together.
    midpoint = int((row.start + row.end) / 2)
    interval_start = int(np.floor(midpoint / ref_interval) * ref_interval)
    interval_end = int(np.ceil(midpoint / ref_interval) * ref_interval)
    gene = f"{row.chr}_{interval_start}_{interval_end}"
    return gene


def cluster_dataframe(df, ref_interval):
    """Process records from tags file."""
    # Create column to keep track of non-assigned genes
    df['no_gene'] = False
    df_no_gene = df.loc[df.gene == '-']
    if len(df_no_gene) > 0:
        # Create a temporary gene name based on chr and location.
        regions = df_no_gene.apply(
            create_region_name, axis=1, args=(ref_interval,))
        df.loc[regions.index, 'gene'] = regions
        df.loc[df.index.isin(regions.index), 'no_gene'] = True
    # Create gene/cell index for subsetting reads prior to clustering.
    df["gene_cell"] = df["gene"] + ":" + df["CB"]
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)
    # UB: corrected UMI tag
    groups = df.groupby("gene_cell")["UR"]
    df["UB"] = groups.transform(cluster)
    df.set_index('read_id', drop=True, inplace=True)
    # Reset unassigned genes to '-'
    df.loc[df.no_gene, 'gene'] = '-'
    df.drop(columns='no_gene', inplace=True)
    return df


class ExpressionSummary(StatsSummary):
    """Gene and transcript feature summary statistics."""

    fields = {
        "tagged",
        "gene_tagged", "transcript_tagged",
        "unique_genes", "unique_transcripts"}

    @classmethod
    def from_pandas(cls, df):
        """Create statistics from pandas dataframe."""
        stats = dict()
        stats["tagged"] = len(df)
        stats["gene_tagged"] = len(df[df.gene != '-'])
        stats["transcript_tagged"] = len(df[df.transcript != '-'])
        stats["genes"] = df['gene'].nunique()  # this will include "-"
        stats["transcripts"] = df['transcript'].nunique()  # and this
        return cls(stats)


def main(args):
    """Run entry point."""
    logger = get_named_logger('CrteMatrix')

    if args.tsv_out is None and args.hdf_out is None:
        raise ValueError("Please supply at least one of `--tsv_out` or `--hdf_out`.")

    logger.info("Reading barcode tag information.")
    df_tags = pd.read_csv(args.barcode_tags, sep='\t', index_col='read_id')

    dups = df_tags[df_tags.index.duplicated(keep='first')]
    if not dups.empty:
        raise ValueError(
            f"One or more input reads are duplicated, please rectify.\n"
            f"Duplicated reads: {list(set(dups.index))[:20]}")

    logger.info("Reading feature information.")
    df_features = pd.read_csv(
        args.features, sep='\t', index_col=0)

    logger.info("Merging barcode and feature information.")
    df_tag_feature = df_tags.merge(
        df_features, how='left', left_index=True, right_index=True).fillna('-')

    logger.info("Filtering reads.")
    df_tag_feature = df_tag_feature.loc[
        (df_tag_feature.CB != '-') & (df_tag_feature.UR != '-')]

    if args.stats:
        logger.info("Writing JSON summary to {args.summary}")
        summary = ExpressionSummary.from_pandas(df_tag_feature)
        summary.to_json(args.stats)

    logger.info("Clustering UMIs.")
    df_tag_feature = cluster_dataframe(df_tag_feature, args.ref_interval)

    logger.info("Preparing output.")
    cols = ['CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'gene', 'transcript', 'start', 'end']
    if len(df_tag_feature) > 0:
        df_tags_out = df_tag_feature[cols].assign(chr=args.chrom)
    else:
        # TODO: comes back to this, it smells janky
        df_tags_out = (
            pd.DataFrame(columns=['read_id'] + cols + ['chr'])
            .set_index('read_id', drop=True)
        )
    df_tags_out.rename(
        columns={v: k for k, v in BAM_TAGS.items()}, copy=False, inplace=True)

    if args.tsv_out:
        logger.info("Writing text output.")
        df_tags_out.to_csv(args.tsv_out, sep='\t')

    if args.hdf_out:
        for feature in ("gene", "transcript"):
            logger.info(f"Creating {feature} expression matrix.")
            matrix = ExpressionMatrix.from_tags(df_tags_out, feature)
            fname = args.hdf_out.with_suffix(f".{feature}{args.hdf_out.suffix}")
            matrix.to_hdf(fname)
