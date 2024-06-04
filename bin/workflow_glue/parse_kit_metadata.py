#!/usr/bin/env python
"""Read in and validate user sample data."""
from pathlib import Path

import pandas as pd

from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("read_samples")
    parent_parser = wf_parser("read_samples_parent")

    parent_parser.add_argument(
        "--kit_config",
        help="Kit-specific details CSV",
        type=Path,
        required=True
    )
    parent_parser.add_argument(
        "--sample_ids",
        help="File with IDs from each sample",
        type=Path,
        required=True
    )
    parent_parser.add_argument(
        "--output",
        help="Output path for merged config",
        type=Path,
        required=True
    )

    subparsers = parser.add_subparsers(help='commands', dest="cmd")

    parser_sheet = subparsers.add_parser(
        'from_sheet', help='Get kit metadata per sample from sample sheet',
        parents=[parent_parser]
    )
    parser_sheet.add_argument(
        "--user_config",
        help="User sample metadata CSV file",
        type=Path,
        required=True
    )

    parser_cli = subparsers.add_parser(
        "from_cli",
        help='Apply the same kit metadata to all samples from CLI variables',
        parents=[parent_parser],
    )
    parser_cli.add_argument(
        "--kit",
        help="10x kit (name:version)",
        required=True
    )
    parser_cli.add_argument(
        "--expected_cells",
        help="Number of expected cells",
        required=True
    )

    return parser


def main(args):
    """Entry point."""
    # Single cell sample sheet expected header
    sc_sample_sheet_header = [
        'sample_id',
        'kit',
        'expected_cells'
    ]

    sample_ids = pd.read_csv(args.sample_ids, index_col=None, header=None)[0].to_list()

    if args.cmd == 'from_cli':
        # No per-sample single-cell sample sheet given by user, so we will use the
        # individual CLI parameters to build a CSV and apply the same parameters to
        # each sample
        entries = [
            [sid.strip(), args.kit, args.expected_cells]
            for sid in sample_ids
        ]
        user_df = pd.DataFrame.from_records(
            entries, columns=sc_sample_sheet_header
        )
    elif args.cmd == 'from_sheet':
        user_df = pd.read_csv(args.user_config)

        # Validate sample sheet header
        if len(set(sc_sample_sheet_header).difference(set(user_df.columns))) != 0:
            raise ValueError(
                'single_cell_sample_sheet should have the following column names: '
                f'{sc_sample_sheet_header}')

    # Validate kit + version combinations.
    kit_df = pd.read_csv(args.kit_config)

    # Check if all supplied kits + version strings are supported
    kit_and_version_diff = set(user_df.kit).difference(kit_df.kit)
    if len(kit_and_version_diff) != 0:
        raise ValueError(
            'the following are not valid kit and version combinations: '
            f'{kit_and_version_diff}')

    # Check that ingressed IDs match sample_ids from sample_sheet
    if set(user_df['sample_id']) != set(sample_ids):
        raise ValueError(
            'Sample IDs from the sc_sample_sheet must match those from those inferred '
            'from the input data:'
            f'\nSample IDs from ingressed data: {sample_ids}'
            f'\nSample IDs from sample sheet: {user_df.sample_id.to_list()}'
            f'\nSamples IDs in sample sheet but not in ingressed data :'
            f'{set(user_df["sample_id"]).difference(sample_ids)}'
            f'\nSamples IDs ingressed data but not in sample sheet :'
            f'{set(sample_ids).difference(user_df["sample_id"])}'
        )

    merged_config = user_df.merge(
        kit_df, on='kit', how='left', suffixes=(None, '_delete'))
    # Create kit name and version columns from the kit:version string
    merged_config[['kit_name', 'kit_version']] \
        = merged_config['kit'].str.split(':', expand=True)
    cols_to_drop = merged_config.columns[merged_config.columns.str.contains('delete')]
    merged_config = merged_config.drop(cols_to_drop, axis=1)
    merged_config.to_csv(args.output, sep=',', index=False)
