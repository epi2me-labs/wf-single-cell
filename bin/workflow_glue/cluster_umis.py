#!/usr/bin/env python
"""Correct UMIs using clustering."""

import itertools
from pathlib import Path

from editdistance import eval as edit_distance
import numpy as np
import pandas as pd
from umi_tools import UMIClusterer

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("cluster_umis")

    parser.add_argument(
        "--chrom",
        help="Chromosome name"
    )

    parser.add_argument(
        "--read_tags",
        help="Read tags TSV file.",
        type=Path
    ),

    parser.add_argument(
        "--output_read_tags",
        help="Output file for read tags TSV. USed for tagging BAM files.",
        type=Path
    )

    parser.add_argument(
        "-i",
        "--ref_interval",
        help="Size of genomic window (bp) to assign as gene name if no gene \
          assigned by featureCounts [1000]",
        type=int,
        default=1000,
    )

    parser.add_argument(
        "--workflow_output",
        help="Output file TSV containing a subset of read-tags in human-readable form",
        type=Path
    )

    parser.add_argument(
        "--feature_assigns",
        help="TSV read gene/transcript assignments file.",
        type=Path
    )

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


def create_region_name(row, ref_interval):
    """
    Create region name.

    Create a 'gene name' based on the aligned chromosome and coordinates.
    The midpoint of the alignment determines which genomic interval to use
    for the 'gene name'.

    :param read: read tags and location
    :type read: tuple
    :param args: object containing all supplied arguments
    :type args: class 'argparse.Namespace'
    :return: Newly created 'gene name' based on aligned chromosome and coords
    :rtype: str
    """
    chrom = row.chr
    start_pos = row.start
    end_pos = row.end

    # Find the midpoint of the alignment
    midpoint = int((start_pos + end_pos) / 2)

    # Pick the genomic interval based on this alignment midpoint
    interval_start = np.floor(midpoint / ref_interval) * ref_interval
    interval_end = np.ceil(midpoint / ref_interval) * ref_interval

    # New 'gene name' will be <chr>_<interval_start>_<interval_end>
    gene = f"{chrom}_{int(interval_start)}_{int(interval_end)}"
    return gene


def cluster(df):
    """Clsuter UMIs.

    Search for UMI clusters within subsets of reads sharing the same corrected barcode
    and gene. In this way the search space is dramatically reduced.

    We are using the UMI-tools directional deduplication method (modified to use
    Levenshtein distance). Connections between nodes within a cluster are generated
    based on edit distance threshold and whether node A counts >= (2* node B counts).
    https://umi-tools.readthedocs.io/en/latest/the_methods.html

    :param df: DataFrame
        Index: gene_cell
        columns: UR (uncorrected UMI), read_id
    :return:
        DataFrame: The same as df with and additional UB (corrected barcode) column.
    """
    # Monkey-patch the umi-tools clusterer with a modified method
    # using Levenshtein instead of Hamming distance.
    UMIClusterer._get_adj_list_directional = get_adj_list_directional_lev
    UMIClusterer.__call__ = umi_clusterer_call

    def umi_tools_cluster(umis):
        clusterer = UMIClusterer(cluster_method="directional")

        umi_map = {}
        umi_counts = umis.value_counts().to_dict()
        clusters = clusterer(umi_counts, threshold=2)

        # Make raw -> corrected umi map
        for clust in clusters:
            if len(clust) == 1:
                # Single UMI cluster. Map it to itself.
                umi_map[clust[0]] = clust[0]
            else:
                for i in range(0, len(clust)):
                    # Map each umi in the cluster to the predicted 'true' umi.
                    umi_map[clust[i]] = clust[0]

        return umis.replace(umi_map)

    # UB: corrected UMI tag
    df["UB"] = df.groupby(
        ["gene_cell"])["UR"].transform(umi_tools_cluster)
    df.set_index('read_id', drop=True, inplace=True)


def process_records(df, args):
    """Process records from tags file.

    For each read, get the gene, barcode and unorrecdted UMI.
    Use that to cluster UMIs to correct errors.
    Write a TSV file including the input column + a corrected UMI tag (UB) column.

    :param df: DataFrame with columns: read_id, UR, gene
    :type df: pd.DataFrame
    """
    df_no_gene = df.loc[df.gene == '-']
    # Create column to keep track of non-assigned genes
    df['no_gene'] = False
    if len(df_no_gene) > 0:
        # For reads with no asssignment, create a temporary gene name
        # based on chr and location. Use that for clustering and reset back to '-' later
        regions = df_no_gene.apply(
            create_region_name, axis=1, args=(args.ref_interval,))
        df.loc[regions.index, 'gene'] = regions
        df.loc[df.index.isin(regions.index), 'no_gene'] = True
    # Create gene/cell index for subsetting reads prior to clustering.
    df["gene_cell"] = df["gene"] + ":" + df["CB"]
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)
    cluster(df)
    # Reset unassigned genes to '-'
    df.loc[df.no_gene, 'gene'] = '-'
    df.drop(columns='no_gene', inplace=True)


def main(args):
    """Run entry point."""
    df_tags = pd.read_csv(args.read_tags, sep='\t', index_col=0)
    df_features = pd.read_csv(
        args.feature_assigns, sep='\t', index_col=0)
    # Merge genes and transcripts onto tags.
    df_tag_feature = df_tags.merge(
        df_features, how='left', left_index=True, right_index=True).fillna('-')

    # Only process reads with a corrected barcode and uncorrected UMI
    df_tag_feature = df_tag_feature.loc[
        (df_tag_feature.CB != '-') & (df_tag_feature.UR != '-')]

    # Process the tag and feature df by adding UB tags inplace.
    process_records(df_tag_feature, args)

    if len(df_tag_feature) > 0:
        # Write a CSV used for tagging BAMs
        df_tags_out = df_tag_feature[
            ['CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'gene', 'transcript']
        ].assign(chr=args.chrom)

        # Write a subset of columns with human-readable names for the final output.
        df_workflow_out = df_tag_feature[
            ['gene', 'transcript', 'CB', 'UB', 'chr', 'start', 'end']]
        df_workflow_out = df_workflow_out.rename(
            columns={
                'CB': 'corrected_barcode',
                'UB': 'corrected_umi'
            }
        )

    else:
        df_tags_out = pd.DataFrame(
            columns=[
                'read_id', 'CR', 'CB', 'CY',
                'UR', 'UB', 'UY', 'gene', 'transcript', 'chr'
             ]
        ).set_index('read_id', drop=True)

        df_workflow_out = pd.DataFrame(
            columns=['read_id', 'gene', 'transcript', 'CB', 'UB', 'chr', 'start', 'end']
        ).set_index('read_id', drop=True)

    df_tags_out.to_csv(args.output_read_tags, sep='\t', index=True)
    df_workflow_out.to_csv(args.workflow_output, sep='\t')
