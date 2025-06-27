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
    parser.add_argument(
        "--ref_interval", type=int, default=1000,
        help="Size of genomic window (bp) to assign as gene name if no gene \
            assigned.")
    parser.add_argument(
        "--chunk_size", type=int, default=200000,
        help="Approximate size of chunks to read in from the input tags file.")
    parser.add_argument(
        "--umi_length", type=int, default=None,
        help="Expected UMI length. Discard reads with corrected UMIs not of this size. \
            If None, no UMI length filtering is done.")
    parser.add_argument(
        "--skip_umi_clustering", action='store_true',
        help="Skip UMI clustering and use pre-corrected UMIs \
             in the UB column of the barcode_tags file.")

    grp = parser.add_argument_group("Output")
    grp.add_argument(
        "--tsv_out", type=Path,
        help="Output TSV containing primary record read-tags in human-readable form")
    grp.add_argument(
        "--sa_tags_out", type=Path,
        help="Output TSV containing supplementary read-tags in human-readable form")
    grp.add_argument(
        "--hdf_out", type=Path,
        help="Output directory for HDF matrix output. \
            HDF chunks will be saved here as <chunk_num>.{gene, transcript}.hdf")
    grp.add_argument(
        "--stats", type=Path,
        help="Output filename for JSON statistics.")

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
    if len(umis) == 1:  # early return; only a single read for this barcode/gene.
        return umis
    clusterer = UMIClusterer(cluster_method="directional")
    umi_counts = collections.Counter(umis)
    clusters = clusterer(umi_counts, threshold=2)

    if len(clusters) == len(umis):  # no corrections, all clusters are singletons
        return umis

    # Create list of corrections
    # The first entry of the cluster is the representative/correct UMI.
    umi_map = {}
    for clust in (x for x in clusters if len(x) > 1):
        correct = clust[0]
        for umi in clust[1:]:
            if umi != correct:
                umi_map[umi] = correct
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


def cluster_dataframe(df, ref_interval, umi_length):
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
    df['gene_cell'] = df['gene'] + ':' + df['CB']
    df['read_id'] = df.index
    df.set_index('gene_cell', inplace=True, drop=True)
    # UB: corrected UMI tag
    groups = df.groupby('gene_cell')['UR']
    df['UB'] = groups.transform(lambda x: cluster(x))
    df.set_index('read_id', drop=True, inplace=True)
    if umi_length is not None:
        valid = df['UB'].str.len() == umi_length
        df.drop(index=df.index[~valid], inplace=True)
    # Reset unassigned genes to '-'
    df.loc[df.no_gene, 'gene'] = '-'
    df.drop(columns='no_gene', inplace=True)
    return df


class ExpressionSummary(StatsSummary):
    """Gene and transcript feature statistics."""

    fields = {
        "valid_barcodes",
        "gene_tagged", "transcript_tagged",
        "unique_genes", "unique_transcripts"}

    def __init__(self, *args, **kwargs):
        """Create summary."""
        super().__init__(*args, **kwargs)
        self.unique_genes = set()
        self.unique_transcripts = set()

    def pandas_update(self, df):
        """Create statistics from pandas dataframe."""
        self["valid_barcodes"] += len(df[df.CB != '-'])
        self["gene_tagged"] += len(df[df.gene != '-'])
        self["transcript_tagged"] += len(df[df.transcript != '-'])
        self.unique_genes.update(df.loc[df.gene != '-', 'gene'])
        self.unique_transcripts.update(df.loc[df.transcript != '-', 'transcript'])

    def to_json(self, fname):
        """Save to JSON. First  convert features counts to n unique."""
        self['genes'] = len(self.unique_genes)
        self['transcripts'] = len(self.unique_transcripts)
        super().to_json(fname)


def chunk_reader(file_path, chunk_size, ub_col=False):
    """
    Read TSV file in chunks, ensuring that no CB (cell barcode) is split between chunks.

    :param file_path str: Path to the TSV
    :param chunk_size int: Desired number of rows per chunk (approximate).
    :param ub_col bool: If True, the 'UB' column is expected to be present.
    :yields: Pandas DataFrame chunk.
    """
    usecols = ['read_id', 'CR', 'CY', 'UR', 'UY', 'chr', 'start', 'end', 'CB', 'SA']
    if ub_col:
        usecols.append('UB')
    reader = pd.read_csv(
        file_path, index_col='read_id', sep="\t", iterator=True,
        usecols=usecols)
    buffer = []
    last_cb = None
    current_cb = None

    # All barcodes must be processed together. A large chunk of records is read in
    # from the barcode pre-sorted DataFrame.
    # We then iterate over rows until a new barcode is encountered,
    # at which point the chunk and the rows are combined and yielded.
    try:
        while True:
            chunk = reader.get_chunk(chunk_size)
            buffer.append(chunk)
            last_cb = chunk.iat[-1, chunk.columns.get_loc("CB")]

            # Look for a new barcode
            while True:
                row = reader.get_chunk(1)
                current_cb = row.at[row.index[0], "CB"]
                if current_cb != last_cb:
                    yield pd.concat(buffer)
                    buffer = [row]  # New barcode for next chunk
                    break
                buffer.append(row)
            last_cb = current_cb
    except StopIteration:
        # No more rows in DataFrame, yield any remaining buffer.
        if buffer:
            yield pd.concat(buffer)


def main(args):
    """Run entry point."""
    logger = get_named_logger('CrteMatrix')

    if args.tsv_out is None and args.hdf_out is None:
        raise ValueError("Please supply at least one of `--tsv_out` or `--hdf_out`.")

    stats = ExpressionSummary()

    tsv_out_cols = [
        'read_id', 'CR', 'CB', 'CY', 'UR', 'UB',
        'UY', 'gene', 'transcript', 'start', 'end', 'chr']
    tag_to_desc_map = {v: k for k, v in BAM_TAGS.items()}
    tsv_out_cols = [tag_to_desc_map.get(t, t) for t in tsv_out_cols]

    # Create header for tags TSV output
    pd.DataFrame(columns=tsv_out_cols).to_csv(
        args.tsv_out, sep='\t', header=True, index=False)

    # Create header for SA tags TSV output
    pd.DataFrame(columns=tsv_out_cols).to_csv(
        args.sa_tags_out, sep='\t', header=True, index=False)

    logger.info("Reading feature information.")
    df_features = pd.read_csv(
        args.features, sep='\t', index_col=0)

    input_has_ub_col = True if args.umi_length is None else False
    for chunk_num, df_tags in enumerate(
            chunk_reader(args.barcode_tags, args.chunk_size, ub_col=input_has_ub_col)):
        logger.info(f'processing chunk: {chunk_num}')

        df_tags = df_tags.merge(
            df_features, how='left', left_index=True, right_index=True).fillna('-')

        logger.info("Filtering reads.")
        df_tags = df_tags.loc[
            (df_tags.CB != '-') & (df_tags.UR != '-')]

        if args.stats:
            logger.info("Writing JSON stats to {args.stats}")
            stats.pandas_update(df_tags)

        if args.skip_umi_clustering:
            logger.info('Skipping UMI clustering.')
        else:
            logger.info("Clustering UMIs.")
            df_tags = cluster_dataframe(df_tags, args.ref_interval, args.umi_length)

        df_tags.rename(
            columns={v: k for k, v in BAM_TAGS.items()}, copy=False, inplace=True)
        df_tags['chr'] = args.chrom
        df_tags.reset_index(inplace=True, drop=False)

        if args.tsv_out:
            logger.info("Writing text output.")
            # Write the tags file for any supplementary records.
            df_sa = df_tags[df_tags.SA != '-']
            df_sa.drop(columns='SA', inplace=True)
            df_sa = df_sa[tsv_out_cols]
            df_sa.to_csv(
                args.sa_tags_out, sep='\t', header=False, mode='a', index=False)
            # Write the tags file for primary records.
            df_tags = df_tags[tsv_out_cols]
            df_tags.to_csv(
                args.tsv_out, sep='\t', header=False, mode='a', index=False)

        if args.hdf_out:
            for feature in ("gene", "transcript"):
                logger.info(f"Creating {feature} expression matrix.")
                matrix = ExpressionMatrix.from_tags(df_tags, feature)
                fname = args.hdf_out / f"{chunk_num}.{feature}.hdf"
                matrix.to_hdf(fname)
    if args.stats:
        stats = stats.to_json(args.stats)
    logger.info("Clustering complete.")
