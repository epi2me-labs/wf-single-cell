import argparse
from itertools import chain
import json
from pathlib import Path
import sys

import pandas as pd
import pytest

import util


ROOT_DIR = Path(__file__).resolve().parent.parent


def args():
    """Parse and process input arguments. Use the workflow params for those missing."""
    # get the path to the workflow output directory
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--input",
        help=(
            "Path to input file / directory with input files / directory with "
            "sub-directories with input files; will take input path from workflow "
            "output if not provided"
        ),
    )
    parser.add_argument(
        "--type",
        choices=util.INPUT_TYPES_EXTENSIONS.keys(),
        help="Input file type",
        required=True,
    )
    parser.add_argument(
        "--wf-output-dir",
        default=ROOT_DIR / "output",
        help=(
            "path to the output directory where the workflow results have been "
            "published; defaults to 'output' in the root directory of the workflow if "
            "not provided"
        ),
    )
    parser.add_argument(
        "--sample_sheet",
        help=(
            "Path to sample sheet CSV file. If not provided, will take sample sheet "
            "path from workflow params (if available)."
        ),
    )
    parser.add_argument(
        "--chunk", type=int,
        help=(
            "Chunk size for output fastq."
        )
    )
    args = parser.parse_args()

    input_type = args.type
    wf_output_dir = Path(args.wf_output_dir)
    ingress_results_dir = (
        wf_output_dir / f"{'xam' if input_type == 'bam' else 'fastq'}_ingress_results"
    )

    # make sure that there are ingress results (i.e. that the workflow has been
    # run successfully and that the correct wf output path was provided)
    if not ingress_results_dir.exists():
        raise ValueError(
            f"{ingress_results_dir} does not exist. Has `wf-template` been run?"
        )

    # get the workflow params
    with open(wf_output_dir / "params.json", "r") as f:
        params = json.load(f)
    input_path = (
        Path(args.input) if args.input is not None else ROOT_DIR / params[input_type]
    )
    sample_sheet = args.sample_sheet
    if sample_sheet is None and params["sample_sheet"] is not None:
        sample_sheet = ROOT_DIR / params["sample_sheet"]

    # Define output type
    output_type = input_type
    if params["wf"]["return_fastq"]:
        output_type = "fastq"

    if not input_path.exists():
        raise ValueError(f"Input path '{input_path}' does not exist.")

    return input_path, input_type, output_type, sample_sheet, ingress_results_dir, args.chunk, params


# prepare data for the tests
@pytest.fixture(scope="module")
def prepare():
    """Prepare data for tests."""
    input_path, input_type, output_type, sample_sheet, ingress_results_dir, chunk_size, params = args()
    valid_inputs = util.get_valid_inputs(input_path, input_type, sample_sheet, chunk_size, params)
    return ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params


# define tests
def test_result_subdirs(prepare):
    """
    Test if workflow results dir contains all expected samples.

    Tests if the published sub-directories in `ingress_results_dir` contain all
    the samples we expect.
    """
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    files = [x for x in ingress_results_dir.iterdir() if x.is_file()]
    subdirs = [x.name for x in ingress_results_dir.iterdir() if x.is_dir()]
    assert not files, "Files found in top-level dir of ingress results"
    assert set(subdirs) == set([meta["alias"] for meta, _ in valid_inputs])


def test_entry_names_and_run_ids(prepare):
    """Test sequence names and run IDs.

    Tests if the concatenated sequences indeed contain all the read IDs of the target
    files in the valid inputs.
    """
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    # unless when run with `--bam ... --wf.return_fastq` the output type (i.e. the type
    # of the files returned by ingress) is the same as the input type
    for meta, path in valid_inputs:
        if path is None:
            # this sample sheet entry had no input dir (or no reads)
            continue
        # get entries in the result file produced by the workflow
        if chunk_size is not None:
            res_seqs_fname = ""
        elif output_type == "fastq":
            res_seqs_fname = "seqs.fastq.gz"
        elif output_type == "bam":
            res_seqs_fname = "reads.bam"
        else:
            raise ValueError(f"Unknown output_type: {output_type}.")

        entries = util.create_preliminary_meta(
            ingress_results_dir / meta["alias"] / res_seqs_fname,
            output_type, chunk_size,
            params["wf"]["return_fastq"])

        # now collect the entries from the individual input files
        exp_read_names = []
        exp_run_ids = []
        target_files = (
            util.get_target_files(path, input_type=input_type)
            if path.is_dir()
            else [path]
        )
        for file in target_files:
            if (
                input_type == "bam"
                and not params["wf"]["keep_unaligned"]
                and util.is_unaligned(file)
            ):
                continue
            curr_entries = util.create_preliminary_meta(
                file, input_type, chunk_size, False)
            exp_read_names += curr_entries["names"]
            exp_run_ids += curr_entries["run_ids"]
        assert set(entries["names"]) == set(exp_read_names)
        assert set(entries["run_ids"]) == set(exp_run_ids)


def test_stats_present(prepare):
    """Tests if the `fastcat` stats are present when they should be."""
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    for meta, path in valid_inputs:
        if path is None:
            # this sample sheet entry had no input dir (or no reads)
            continue
        if output_type == "fastq":
            expect_stats = params["wf"]["fastcat_stats"]
            stats_dir_name = "fastcat_stats"
            stats_file_names = [
                "per-file-stats.tsv",
                "per-read-stats.tsv.gz",
                "run_ids",
                "length.hist",
                "quality.hist"
            ]
        else:
            # `bamstats` we only expect when they were requested
            expect_stats = params["wf"]["bamstats"]
            stats_dir_name = "bamstats_results"
            stats_file_names = [
                "bamstats.readstats.tsv.gz",
                "bamstats.flagstat.tsv",
                "run_ids",
                "accuracy.hist",
                "coverage.hist",
                "length.hist",
                "quality.hist"
            ]
        stats_dir = ingress_results_dir / meta["alias"] / stats_dir_name
        # assert that stats are there when we expect them
        assert expect_stats == stats_dir.exists()
        # make sure that the per-file stats, per-read stats, and run ID files are there
        if expect_stats:
            for fname in stats_file_names:
                assert (
                    ingress_results_dir / meta["alias"] / stats_dir_name / fname
                ).is_file()


def test_metamap(prepare):
    """Test if the metamap in the ingress results is as expected."""
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    for meta, _ in valid_inputs:
        # prepare() uses a function to parse both inputs and outputs,
        # add in some output specific things
        meta = util.add_output_n_fastq(meta, output_type, chunk_size)
        sample_results = ingress_results_dir / meta["alias"]

        # if there were no stats, we can't expect run IDs in the metamap
        if not list(sample_results.glob("*stats*/run_ids")):
            meta["run_ids"] = []
            # if there are not stats, reset extra fields to defaults
            # could not be stats where barcodes are in sample sheet but not in data.
            if output_type == "fastq":
                meta["n_seqs"] = None
            elif output_type == "bam":
                meta["n_primary"] = None
                meta["n_unmapped"] = None

        # read what nextflow had
        with open(sample_results / "metamap.json", "r") as f:
            metamap = json.load(f)
        assert meta == metamap


def test_reads_sorted(prepare):
    """If input type is BAM, test if the emitted files were sorted."""
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    if input_type == "fastq":
        return
    for meta, _ in valid_inputs:
        stats_file = (
            ingress_results_dir
            / meta["alias"]
            / "bamstats_results"
            / "bamstats.readstats.tsv.gz"
        )
        if stats_file.exists():
            stats_df = pd.read_csv(stats_file, sep="\t", index_col=0)
            # check that the start coordinates of all aligned reads are sorted within
            # their respective reference
            assert (
                stats_df.query('ref != "*"')
                .groupby("ref")["rstart"]
                .is_monotonic_increasing.all()
            )


def test_reads_index(prepare):
    """If input type is BAM, check that the BAI index exists."""
    ingress_results_dir, input_type, output_type, valid_inputs, chunk_size, params = prepare
    if output_type == "fastq":
        return
    for meta, path in valid_inputs:
        if path is None:
            # this sample sheet entry had no input dir (or no reads)
            continue
        # Create BAI file path
        bai_file = (
            ingress_results_dir
            / meta["alias"]
            / 'reads.bam.bai'
        )
        if not bai_file.is_file():
            raise ValueError(f"Missing index: {bai_file.as_posix()}.")


if __name__ == "__main__":
    # trigger pytest
    ret_code = pytest.main([Path(__file__).resolve(), "-vv", "-s"])
    sys.exit(ret_code)
