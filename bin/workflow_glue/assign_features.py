"""Assign reference genes and transcripts to read_id."""
from pathlib import Path as Path
import re

import pandas as pd
from pysam import AlignmentFile

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("assign_features")

    parser.add_argument(
        "transcriptome_bam", type=Path,
        help="Bam file from alignment of reads to the assembled transcriptome.")
    parser.add_argument(
        "gffcompare_tmap", type=Path,
        help="The .tmap output from gffcompare")
    parser.add_argument(
        "gtf", type=Path,
        help="GTF for chromosome. Should contain gene_id and gene_name.")
    parser.add_argument(
        "tags", type=Path,
        help="TSV file with read_id index and genomic alignment mapq column.")
    parser.add_argument(
        "output", type=Path,
        help="Path to output file")
    parser.add_argument(
        "--chunksize", type=int, default=100000,
        help="Process the BAM in chunks no larger than this.")
    parser.add_argument(
        "--min_tr_coverage", type=float, default=0.4,
        help="Minimum transcript coverage for assignment.")
    parser.add_argument(
        "--min_read_coverage", type=int, default=0.4,
        help="Minimum read coverage for assignment.")
    parser.add_argument(
        "--min_mapq", type=int, default=30,
        help="mapq threshold from genomic alignment above which reads will be assigned "
             "a gene and transcript")
    return parser


def parse_gtf(annotation_file):
    """Get a mapping of query transcript id to reference gene name."""
    with open(annotation_file, 'r') as fh:

        gene_regex = re.compile(r'(gene_name|genename) "(.+?)"')
        tr_regex = re.compile(r'transcript_id "(.+?)"')
        records = []
        for line in fh:
            if line.startswith('#'):
                continue
            if line.split('\t')[2].lower() not in ['transcript', 'mrna']:
                continue
            g = gene_regex.search(line, re.IGNORECASE)
            if g:
                gene_name = g.group(2)
            else:
                gene_name = None
            t = tr_regex.search(line, re.IGNORECASE)
            if t:
                transcript_id = t.group(1)
            else:
                transcript_id = None

            records.append([transcript_id, gene_name])
        df = pd.DataFrame.from_records(records, columns=['ref_id', 'gene_name'])
        return df


def parse_bam(transcriptome_bam, chunksize):
    """Create batched dataframes summarising alignments.

    Alignments with equal query name are guaranteed to be grouped in one batch
    provided the input file is namesorted.

    :param chunksize: (approximate) number of alignments to yield in one go.
    """
    cols = [
        'read_id', 'query_transcript', 'aln_score', 'tr_cov', 'q_cov', 'tr_mapq']

    with AlignmentFile(transcriptome_bam, "rb", check_sq=False) as bam:
        ref_lengths = dict(zip(bam.header.references, bam.header.lengths))
        records = []
        last_name = None
        for align in bam.fetch(until_eof=True):
            if align.reference_name is None:
                continue
            current_name = align.query_name
            try:
                aln_score = align.get_tag('AS')
            except KeyError:
                aln_score = 0
            map_st = align.reference_start
            map_en = align.reference_end
            tr_cov = float(map_en - map_st) / ref_lengths[align.reference_name]

            tr_mapq = align.mapping_quality
            # Query transcript built by stringtie is the refernce here
            query_transcript = align.reference_name
            query_coverage = float(
                align.query_alignment_length / align.infer_read_length())
            # The chunk size limit has been reached. Make sure alignemnts from
            # the same read are not split into seperate chunks
            if len(records) >= chunksize and current_name != last_name:  # doc
                df = pd.DataFrame.from_records(records, columns=cols)
                records = []
                yield df

            records.append([
                align.query_name, query_transcript,
                aln_score, tr_cov, query_coverage, tr_mapq
            ])

            last_name = align.query_name

    if len(records) > 0:
        df = pd.DataFrame.from_records(records, columns=cols)
        yield df


def main(args):
    """Assign gene and transcript to reads.

    Parse the transcriptome-aligned BAM and attempt to disambiquate reads that map
    to multiple transcripts.

    Filtering criteria borrowed from FLAMES:
    https://github.com/LuyiTian/FLAMES/blob/774e16ae53a1430e03081970827e93f1fbaecead/
    python/count_tr.py#L101
    """
    logger = get_named_logger('AssignFeat')
    logger.info('Assigning genes and transcripts to reads.')

    logger.info("Loading mapping quality from genomic alignment.")
    df_genomic_mapq = pd.read_csv(
        args.tags, sep='\t', index_col=0, usecols=['read_id', 'mapq'])
    df_genomic_mapq.rename(columns={'mapq': 'genome_mapq'}, inplace=True)

    # Load gffcompare output that maps query transcript to reference transcript
    logger.info("Loading gffcompare output.")
    df_gffcompare_tmap = pd.read_csv(
        args.gffcompare_tmap, index_col=None, sep='\t')

    logger.info("Loading annotation.")
    df_ann = parse_gtf(args.gtf)

    # Write dataframe header to file
    pd.DataFrame(
        columns=['gene', 'transcript']).to_csv(
        args.output, sep='\t', header=True)

    logger.info("Processing BAM.")
    total_alignments = 0
    for df in parse_bam(args.transcriptome_bam, args.chunksize):
        total_alignments += len(df)

        # Merge the read alignments with the gffcompare tmap output to
        # assign read_id to the reference transcript.
        # Any entries in the tmap file that do not have a ref_id (reference transcript)
        # are potentially novel but are currently unhandled.
        df_tr = df.merge(
            df_gffcompare_tmap[['qry_id', 'ref_id', 'class_code']],
            left_on='query_transcript',
            right_on='qry_id')

        # Merge the gene names from the original input annotation gtf.
        # This merging discards reads that were not present in the input tags file
        # These would have been discarded at earlier parts of the workflow.
        df_tr = df_tr.merge(df_ann, how='inner', left_on='ref_id', right_on='ref_id')
        df_tr.gene_name = df_tr.gene_name.fillna('-')
        df_tr.rename(
            columns={
                'ref_id': 'transcript',
                'gene_name': 'gene',
            }, inplace=True)

        # Remove unwanted categories of transcripts:
        # https://ccb.jhu.edu/software/stringtie/gffcompare.shtml
        df_tr.loc[df_tr['class_code'].isin(['i', 'y', 'p', 's']), 'transcript'] = '-'

        df_tr = df_tr.merge(
            df_genomic_mapq, how='left', left_on='read_id', right_on='read_id')

        multimaped_idxs = df_tr.read_id.duplicated(keep=False)
        df_multimap = df_tr.loc[multimaped_idxs]
        df_uniqmap = df_tr.loc[~multimaped_idxs]

        # Process the uniquely-mapping isoforms
        df_uniqmap.loc[df_uniqmap.genome_mapq < args.min_mapq, 'gene'] = '-'
        df_uniqmap = df_uniqmap.loc[
            df_uniqmap.tr_mapq > 0][['read_id', 'gene', 'transcript']]

        # Proces the multi-mapping reads
        # find the best transcript per read.
        # pseudo code:
        #   for each read:
        #       # to assign a gene
        #       sort aligns by mapq and take best
        #       # to assign a transcript
        #       if (best aln score unique) or (q_cov of best aln unique):
        #           and if tr_cov of best aln > min_tr_coverage
        #           and if q_cov of best aln > min_read_coverage
        #               then assign best aln transcript
        #       elif (tr_cov of best aln > 0.8 and tr_cov of best aln unique):
        #           assign best aln transcript
        #       else
        #           assign "-"

        df_multimap = df_multimap.sort_values(
            ['read_id', 'aln_score', 'q_cov', 'tr_cov'], ascending=False)

        # find the best gene per read. Note this corrects two bugs in the previous code:
        #    1. previous was taking worst mapq: .sort_value().iloc[0]
        #    2. previous was score > minmapq, but should be >=
        read_groups = df_multimap.groupby('read_id', sort=False)
        idx = read_groups['genome_mapq'].idxmax()
        best_gene = df_multimap.loc[idx].set_index('read_id')
        best_gene.loc[best_gene['genome_mapq'] < args.min_mapq, "gene"] = "-"
        best_gene = best_gene[["gene"]]

        # calculate diffs between consecutive rows, we only care about the
        # first difference
        df_multimap['as_diff'] = read_groups["aln_score"].diff() != 0
        df_multimap['q_cov_diff'] = read_groups["q_cov"].diff() != 0
        df_multimap['tr_cov_diff'] = read_groups["tr_cov"].diff() != 0
        stats = read_groups.nth(1).set_index('read_id')[
            ['as_diff', 'q_cov_diff', 'tr_cov_diff']]

        # calculate conditions for the if
        stats["cond1"] = stats['as_diff'] | stats['q_cov_diff']

        # the and conditions for the if
        stats["good_tr_cov"] = (
            read_groups.nth(0).set_index('read_id')["tr_cov"] >= args.min_tr_coverage)
        stats["good_read_cov"] = (
            read_groups.nth(0).set_index('read_id')["q_cov"] >= args.min_read_coverage)
        stats["cond2"] = stats["good_tr_cov"] & stats["good_read_cov"]

        # conditions for the elif
        stats["ok_tr_cov"] = read_groups.nth(0).set_index('read_id')["tr_cov"] > 0.8
        stats["cond3"] = stats['tr_cov_diff'] & stats['ok_tr_cov']

        # combine the if-elif conditions
        keep = (stats["cond1"] & stats["cond2"]) | (~stats["cond1"] & stats["cond3"])
        stats["transcript"] = read_groups.nth(0).set_index('read_id')["transcript"]
        stats.loc[~keep, "transcript"] = "-"
        stats = stats[["transcript"]]

        # create a final table with the best gene and transcript
        df_demultimapped = stats.join(best_gene).reset_index(drop=False)

        # tidy up and write out the chunk
        df_assigned = pd.concat([df_uniqmap, df_demultimapped])
        df_assigned = df_assigned.loc[
            (df_assigned.gene != '-') | (df_assigned.transcript != '-')]
        df_assigned.set_index('read_id', drop=True, inplace=True)
        df_assigned = df_assigned.fillna('-')[['gene', 'transcript']]
        df_assigned.to_csv(args.output, mode='a', sep='\t', header=False)
        logger.info(f"Processed {total_alignments} alignments.")
    logger.info("Finished.")
