#!/usr/bin/env python
"""Correct UMIs using clustering."""

import itertools
from pathlib import Path

from editdistance import eval as edit_distance
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
        "--workflow_output",
        help="Output file TSV containing a subset of read-tags in human-readable form",
        type=Path
    )

    parser.add_argument(
        "--cell_gene_max_reads",
        help="Maximum number of reads to consider for a particular \
        gene + cell barcode combination. \
        This is required to prevent too many PCR \
        duplicates from crashing the UMI clustering algorithm. \
        Can be increased \
        if sufficient UMI complexity is observed. [20000]",
        type=int,
        default=20000
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
    return df.set_index('read_id', drop=True)


def process_records(df):
    """Process records from tags file.

    For each read, get the gene, barcode and unorrecdted UMI.
    Use that to cluster UMIs to correct errors.
    Write a TSV file including the input column + a corrected UMI tag (UB) column.

    :param df: DataFrame with columns: read_id, UR, gene
    :type df: pd.DataFrame
    """
    # Only process recoreds with gene, corrected barcode and uncorrected UMIs.
    df = df.loc[(df.gene != '-') & (df.CB != '-') & (df.UR != '-')]

    # Create gene/cell index for subsetting reads prior to clustering.
    df["gene_cell"] = df["gene"] + ":" + df["CB"]
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)

    df_ub = cluster(df)
    return df_ub


def main(args):
    """Run entry point."""
    df_tags = pd.read_csv(args.read_tags, sep='\t', index_col=0)
    df_features = pd.read_csv(
        args.feature_assigns, sep='\t', index_col=0)  # ?keep_default_na=False
    # Merge genes and transcripts onto tags.
    df_tag_feature = df_features.merge(df_tags, left_index=True, right_index=True)

    df_ub = process_records(df_tag_feature)
    if len(df_ub) > 0:
        # Merge The corrected UMIs back onto the orignal dataframe
        # What happend to non-assigned rows?
        df_ub = df_tag_feature.merge(
            df_ub[['UB']],
            how='left', left_index=True, right_index=True)
        df_ub.UB.fillna('-', inplace=True)

        # Write a CSV used for tagging BAMs
        df_tags_out = df_ub[
            ['CR', 'CB', 'CY', 'UR', 'UB', 'UY', 'gene', 'transcript']
        ].assign(chr=args.chrom)

        # Write a subset of columns with human-readable names for the final output.
        df_workflow_out = df_ub[
            ['gene', 'transcript', 'CB', 'UB', 'chr', 'start', 'end']]
        df_workflow_out.rename(
            columns={
                'CB': 'corrected_barcode',
                'UB': 'corrected_umi'
            }, inplace=True
        )

    else:
        df_tags_out = pd.DataFrame()
        df_workflow_out = pd.DataFrame()

    df_tags_out.to_csv(args.output_read_tags, sep='\t', index=True)
    df_workflow_out.to_csv(args.workflow_output, sep='\t')
