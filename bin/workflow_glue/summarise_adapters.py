#!/usr/bin/env python
"""Summarise the adapter configurations."""
from collections import Counter
import json
from pathlib import Path

import numpy as np
import pandas as pd
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("summarise_adapters")

    parser.add_argument(
        "--read_config_tsv",
        help="Path to file with adapter configs",
        type=Path
    )
    parser.add_argument(
        "--out",
        help="Path to output JSON file",
        type=Path
    )
    parser.add_argument(
        "--sample_id",
        help="ID of the sample"
    )
    parser.add_argument(
        "--threads",
        help="Number of threads to use",
        type=int,
        default=2
    )
    return parser


def main(args):
    """Run entry point.

    Make a JSON file with read summary statistics and counts of the various adapter
    configurations.
    """
    logger = get_named_logger('SummAdapts')
    logger.info("Summarising adapter configurations")

    configs = Counter()
    labels = Counter()
    n_full_lengths = 0
    n_stranded = 0
    n_plus = 0
    n_minus = 0
    chunk_sizes = []
    mean_lengths = []
    variances = []

    with pd.read_csv(
            args.read_config_tsv, chunksize=50000,
            sep='\t',
            usecols=[
                'read_id', 'readlen', 'fl', 'lab', 'orig_strand',
                'orig_adapter_config']) as reader:

        for df in reader:
            chunk_sizes.append(len(df))
            mean_lengths.append(df['readlen'].sum())
            variances.append(np.var(df['readlen']))
            n_full_lengths += len(df[df.fl].index)
            n_stranded += len(df[df['orig_strand'] != '*'])
            n_plus += len(df[df['orig_strand'] == '+'])
            n_minus += len(df[df['orig_strand'] == '-'])
            configs.update(dict(df['orig_adapter_config'].value_counts()))
            labels.update(dict(df['lab'].value_counts()))

    # Calculate pooled variance and standard deviation from chunks
    c = np.array(chunk_sizes)
    v = np.array(variances)
    var = sum((c - 1) * v) / sum(c - 1)
    stdev = np.sqrt(var)

    stats = {}
    stats[args.sample_id] = {}
    stats[args.sample_id]["general"] = {}
    stats[args.sample_id]["general"]["n_reads"] = sum(chunk_sizes)
    # Chunks should be similar size so no need for a weighted average
    stats[args.sample_id]["general"]["rl_mean"] = (
            np.sum(mean_lengths) / sum(chunk_sizes))
    stats[args.sample_id]["general"]["rl_std_dev"] = stdev
    stats[args.sample_id]["general"]["n_fl"] = n_full_lengths
    stats[args.sample_id]["general"]["n_stranded"] = n_stranded

    stats[args.sample_id]["strand_counts"] = {}
    stats[args.sample_id]["strand_counts"]["n_plus"] = n_plus
    stats[args.sample_id]["strand_counts"]["n_minus"] = n_minus
    stats[args.sample_id]["detailed_config"] = {}

    for category, n in configs.items():
        stats[args.sample_id]["detailed_config"][category] = int(n)

    stats[args.sample_id]["summary_config"] = {}
    for label, n in labels.items():
        stats[args.sample_id]["summary_config"][label] = int(n)

    with open(args.out, "w") as f:
        json.dump(stats, f, indent=4)
