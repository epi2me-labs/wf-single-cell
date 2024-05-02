"""Aggregate adapter configuration summaries."""
from pathlib import Path

from .adapter_scan_vsearch import AdapterSummary  # noqa: ABS101
from .util import get_named_logger, wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("summarise_adapters")

    parser.add_argument(
        "input_dir", type=Path,
        help="Path to JSON files to aggregate.")
    parser.add_argument(
        "output", type=Path,
        help="Path to output JSON file")

    return parser


def main(args):
    """Aggregate multiple adapter configuration summary files."""
    logger = get_named_logger('AggAdptCnf')
    logger.info("Aggregating adapter configurations")

    fnames = list(args.input_dir.glob("*.json"))
    if len(fnames) == 0:
        raise IOError("No summary JSON files found.")

    summary = AdapterSummary.from_json(fnames[0])
    if len(fnames) > 1:
        for other in fnames[1:]:
            summary += AdapterSummary.from_json(other)
    summary.to_json(args.output)
