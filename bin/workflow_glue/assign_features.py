#!/usr/bin/env python
"""Assign reference genes and transcripts to read_id."""
import re

import pandas as pd
from pysam import AlignmentFile

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("assign_features")

    parser.add_argument(
        "--transcriptome_bam",
        help="Bam file from alignment of reads to the assembled transcriptome."
    )

    parser.add_argument(
        "--gffcompare_tmap",
        help="The .tmap output from gffcompare"
    )

    parser.add_argument(
        "--gtf",
        help="GTF for chromosome. Should contain gene_id and gene_name "
    )

    parser.add_argument(
        "--tags",
        help="TSV file with read_id index and genomic alignment mapq column."
    )

    parser.add_argument(
        "--output",
        help="Path to output file"
    )

    parser.add_argument(
        "--chunksize",
        type=int,
        default=10000,
        help="Process the BAM in chunks no larger than this."
    )

    parser.add_argument(
        "--min_tr_coverage",
        type=int,
        default=0.4,
        help="Minimum transcript coverage for assignment."
    )

    parser.add_argument(
        "--min_read_coverage",
        type=int,
        default=0.4,
        help="Minimum read coverage for assignment."
    )

    parser.add_argument(
        "--min_mapq",
        help="mapq threshold from genomic alignment above which reads will be assigned "
             "a gene and transcript",
        type=int, default=30
    )
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
    """Parse the transcriptome alignment BAM sorted by name.

    Yield chunks with a maximum size 100,000
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

    # Load genomic alignment mapq scores for filtering gene calls
    df_genomic_mapq = pd.read_csv(
        args.tags, sep='\t', index_col=0, usecols=['read_id', 'mapq'])
    df_genomic_mapq.rename(columns={'mapq': 'genome_mapq'}, inplace=True)

    # Load gffcompare output that maps query transcript to reference transcript
    df_gffcompare_tmap = pd.read_csv(
        args.gffcompare_tmap,
        index_col=None,
        sep='\t')

    df_ann = parse_gtf(args.gtf)

    # Write dataframe header to file
    pd.DataFrame(
        columns=['gene', 'transcript']).to_csv(
        args.output, sep='\t', header=True)

    for df in parse_bam(args.transcriptome_bam, args.chunksize):

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

        assigned = []

        # Attempt to choose one of the multi-mapped alignments
        for read_id, df_read in df_multimap.groupby('read_id'):
            df_read.sort_values(
                ['aln_score', 'q_cov', 'tr_cov'], ascending=False, inplace=True)

            tr = '-'
            # Choose alignment with the best alignment score or query coverage
            if (
                    df_read.iloc[0].aln_score > df_read.iloc[1].aln_score
                    or df_read.iloc[0].q_cov > df_read.iloc[1].q_cov
            ):
                # Assign if there is enough read coverage
                if df_read.iloc[0].tr_cov >= args.min_tr_coverage and \
                        df_read.iloc[0].q_cov >= args.min_read_coverage:
                    tr = df_read.iloc[0].transcript

            # If there's an alignment score (AS) and read coverage tie,
            # but a higher transcript coverage
            elif df_read.iloc[0].tr_cov > df_read.iloc[1].tr_cov:
                if df_read.iloc[0].tr_cov > 0.8:
                    tr = df_read.iloc[0].transcript

            # Assign gene
            top_gene = df_read.sort_values('genome_mapq').iloc[0]
            if top_gene.genome_mapq > args.min_mapq:
                gene = top_gene.gene
            else:
                gene = '-'
            assigned.append([read_id, tr, gene])

        # Merge the resolved multimapped alignments back with the unique mapper
        df_demultimapped = pd.DataFrame.from_records(
            assigned, columns=['read_id', 'transcript', 'gene'])
        df_assigned = pd.concat([df_uniqmap, df_demultimapped])

        df_assigned = df_assigned.loc[
            (df_assigned.gene != '-') | (df_assigned.transcript != '-')]

        df_assigned.set_index('read_id', drop=True, inplace=True)

        df_assigned = df_assigned.fillna('-')[['gene', 'transcript']]

        # Write the processed chunk
        df_assigned.to_csv(args.output, mode='a', sep='\t', header=False)
