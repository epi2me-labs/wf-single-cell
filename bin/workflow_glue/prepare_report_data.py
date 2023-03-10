#!/usr/bin/env python
"""Preapre data for the report."""
from collections import defaultdict
import json

import pandas as pd

from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("prepare_report_data")

    parser.add_argument(
        "--sample_id",
        help="ID of the sample being processed")
    parser.add_argument(
        "--read_tags",
        help="TSV with read id gene, transcript, barcode, umi assignments")
    parser.add_argument(
        "--config_stats",
        help="Workflow sumamry statistics")
    parser.add_argument(
        "--white_list",
        help="Workflow sumamry statistics")
    parser.add_argument(
        "--output", help="Output csv with workflow data.")
    return parser


def _get_sample_summaries(read_tags, white_list):

    df = pd.read_csv(read_tags, sep='\t')
    total_cells = len(pd.read_csv(white_list, sep='\t'))
    # The df contains only reads with a barcode and umi tag
    total_tagged = len(df)
    gene_tagged_df = df[df.gene != '-']
    gene_tagged = len(gene_tagged_df)
    total_genes = len(gene_tagged_df['gene'].unique())
    transcript_tagged_df = df[df.transcript != '-']
    transcript_tagged = len(transcript_tagged_df)
    total_transcripts = len(transcript_tagged_df['transcript'].unique())
    cols = [
        'total_tagged', 'gene_tagged', 'transcript_tagged',
        'total_genes', 'total_transcripts', 'total_cells']
    record = [[
        total_tagged, gene_tagged, transcript_tagged,
        total_genes, total_transcripts, total_cells]]

    df = pd.DataFrame.from_records(record, columns=cols)
    return df


def _parse_config_stats(wf_config_stats, summ_df, sample_id):
    """Parse the workflow config_stats."""
    results = defaultdict(dict)
    with open(wf_config_stats) as json_file:
        data = json.load(json_file)
        gen = data[sample_id]['general']
        summ = data[sample_id]['summary_config']

        # general section of the json
        general_keys = [
            'n_reads',
            'n_fl',
            'n_stranded',
            'n_plus',
            'n_minus']
        for g in general_keys:
            val = gen.get(g, 0)
            results[sample_id][g] = val

        # Summary section
        summ_keys = [
            'double_adapter2',
            'single_adapter2',
            'single_adapter1',
            'double_adapter1',
            'no_adapters',
            'other']
        for s in summ_keys:
            val = summ.get(s, 0)
            results[sample_id][s] = val
        df = pd.DataFrame.from_dict(results).reset_index()
        s = summ_df.T.reset_index(drop=False)
        s.columns = ['index', sample_id]
        df = pd.concat([df, s])

        df = df.melt(id_vars='index')
        # Class can be workflow stage or the primer config.
        df = df.rename(
            columns={
                'index': 'class',
                'variable': 'sample',
                'value': 'count'})  # Proportion of reads [%]
        df['class'].replace('n_fl', 'full_length', inplace=True)
        dfs = []
        for _, df_sample in df.groupby('sample'):
            df_sample[
                'Proportion of reads [%]'] = 100 / df_sample.loc[
                    df_sample['class'] == 'n_reads', 'count'].values[0] \
                    * df_sample['count']
            dfs.append(df_sample)
        final_df = pd.concat(dfs)
        return final_df


def main(args):
    """Entry point for script."""
    logger = get_named_logger('PrepReport')
    logger.info('preparing report data')
    df_summ = _get_sample_summaries(
        args.read_tags, args.white_list)

    df_survival = _parse_config_stats(
        args.config_stats, df_summ, args.sample_id)

    df_survival['sample_id'] = args.sample_id
    df_summ['sample_id'] = args.sample_id

    # Put n_reads into the summary df
    n_reads = df_survival.loc[
        df_survival['class'] == 'n_reads', 'count'].values[0]
    df_summ['n_reads'] = n_reads
    df_summ.to_csv('sample_summary.tsv', sep='\t')
    df_survival.to_csv('survival_data.tsv', sep='\t')


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
