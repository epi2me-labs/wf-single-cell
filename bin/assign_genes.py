import argparse
import collections
import logging
import os
import pathlib
import shutil
import sys

import bioframe as bf
import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)


def parse_args():
    # Create argument parser
    parser = argparse.ArgumentParser()

    # Positional mandatory arguments
    parser.add_argument(
        "bed",
        help="BED file of alignments intervals",
        type=str,
    )

    parser.add_argument(
        "gtf",
        help="GTF file of gene annotations",
        type=str,
    )

    # Optional arguments
    parser.add_argument(
        "-q",
        "--mapq",
        help="Minimum mapping quality to use for feature assignment [60]",
        type=int,
        default=60,
    )

    parser.add_argument(
        "--output",
        help="Output file [./read_annotations.tsv]",
        type=str,
        default="./read_annotations.tsv",
    )

    parser.add_argument(
        "-c",
        "--chunk_size",
        help="BED alignments per chunk to process [200000]",
        type=int,
        default=200000,
    )

    parser.add_argument(
        "--verbosity",
        help="logging level: <=2 logs info, <=3 logs warnings",
        type=int,
        default=2,
    )

    # Parse arguments
    args = parser.parse_args()

    return args


def init_logger(args):
    logging.basicConfig(
        format="%(asctime)s -- %(message)s", datefmt="%Y-%m-%d %H:%M:%S"
    )
    logging_level = args.verbosity * 10
    logging.root.setLevel(logging_level)
    logging.root.handlers[0].addFilter(lambda x: "NumExpr" not in x.msg)


def load_gtf(args):
    """
    Load the GTF file of reference gene annotations.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: dataframe of gene annotation intervals
    :rtype: pandas dataFrame
    """
    cols = [
        "chrom",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    df = pd.read_csv(args.gtf, sep="\t", comment="#", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "GTF file not loading as a valid dataframe!"

        # Restrict annotations to the gene level and isolate gene name from attributes
        df = df[df["feature"] == "gene"]
        df["attribute"] = df["attribute"].str.split(";", expand=True).iloc[:, 3]
        df["attribute"] = df["attribute"].str.split(" ", expand=True).iloc[:, -1]
        df["attribute"] = df["attribute"].str.replace('"', "")

    return df


def load_bed(args):
    """
    Read in the BED file of alignment intervals. BED file created by running
    bedtools bamtobed -i <BAM>

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: dataframe of alignment intervals
    :rtype: pandas dataFrame
    """
    cols = [
        "chrom",
        "start",
        "end",
        "name",
        "score",
        "strand",
    ]
    df = pd.read_csv(args.bed, sep="\t", header=None, names=cols)
    if df.shape[0] > 0:
        assert bf.is_bedframe(df), "BED file not loading as a valid dataframe!"

    df["aln_len"] = df["end"] - df["start"]

    return df


def assign_status_low_mapq(df, args):
    """
    Assign all reads with ambiguous mapping as "Unassigned_mapq" and set the
    assigned gene to NA.

    :param df: dataframe of all gene annotation assignments per alignment
    :type df: pandas dataFrame
    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: dataframe of gene annotation assignments with Unassigned_mapq labels
    :rtype: pandas dataFrame
    """
    df.loc[df["score"] < args.mapq, "status"] = "Unassigned_mapq"
    df.loc[df["score"] < args.mapq, "gene"] = "NA"
    return df


def assign_status_ambiguous_overlap(df):
    """
    Assign all reads with equal overlap with multiple features as
    "Unassigned_ambiguous" and set gene assignment to NA

    :param df: dataframe of all gene annotation assignments per alignment
    :type df: pandas dataFrame
    :return: dataframe of gene annotation assignments with Unassigned_ambiguous labels
    :rtype: pandas dataFrame
    """
    is_ambiguous = df[df["status"] == "Unknown"].duplicated(
        subset=["index_bed", "overlap_bp"], keep=False
    )
    ambiguous_idx = is_ambiguous.index[is_ambiguous]
    df.loc[ambiguous_idx, "status"] = "Unassigned_ambiguous"
    df.loc[ambiguous_idx, "gene"] = "NA"
    df = df.drop_duplicates(subset=["index_bed", "overlap_bp", "status"])

    return df


def assign_status_no_features(df):
    """
    Assign all reads without any features as "Unassigned_no_features"

    :param df: dataframe of all gene annotation assignments per alignment
    :type df: pandas dataFrame
    :return: dataframe of gene annotation assignments with Unassigned_no_features labels
    :rtype: pandas dataFrame
    """
    df.loc[df["gene"] == 0, "status"] = "Unassigned_no_features"
    df.loc[df["gene"] == 0, "gene"] = "NA"
    return df


def find_largest_overlap(df):
    """
    Find the largest overlap when multiple overlaps are present. Label the longest
    overlap as Assigned and drop the shorter overlaps from the dataframe.

    :param df: dataframe of all gene annotation assignments per alignment
    :type df: pandas dataFrame
    :return: dataframe of gene annotation assignments with Assigned label
    :rtype: pandas dataFrame
    """
    # Find the indices of the largest overlap for each alignment entry
    max_ovlp_idx = df.groupby(["index_bed"])["overlap_bp"].idxmax().sort_values().values

    # Only keep the entry with the largest overlap
    df = df.loc[max_ovlp_idx, :]
    unknown_idx = [i for i in max_ovlp_idx if df.loc[i, "status"] == "Unknown"]
    df.loc[unknown_idx, "status"] = "Assigned"
    return df


def get_overlaps(bed, gtf):
    """
    Use bioframe to identify overlaps between the alignment and gene annotation
    intervals.

    :param bed: dataframe containing alignment intervals
    :type bed: pandas dataFrame
    :param gtf: dataframe containing gene annotation intervals
    :type gtf: pandas dataFrame
    :return: dataframe of annotations per alignment based on overlaps
    :rtype: pandas dataFrame
    """

    # how="left" means that the returned dataframe contains all alignments,
    # regardless of whether it overlapped with any annotation intervals.
    df = bf.overlap(
        bed,
        gtf,
        how="left",
        suffixes=("_bed", "_gtf"),
        return_overlap=True,
        return_index=True,
    )

    # Now keeping index_bed to keep track of the input order. We'll want to
    # get the largest overlap for each index_bed value, which will properly
    # handle the issue of supplementary alignments in the BAM.

    # Only keep relevant columns
    df = df[
        [
            "index_bed",
            "name_bed",
            "chrom_bed",
            "score_bed",
            "strand_gtf",
            "attribute_gtf",
            "overlap_start",
            "overlap_end",
        ]
    ].fillna(0)

    df = df.rename(
        columns={
            "name_bed": "read",
            "chrom_bed": "chrom",
            "score_bed": "score",
            "strand_gtf": "strand",
            "attribute_gtf": "gene",
        }
    )
    df["score"] = df["score"].astype(int)
    df["overlap_bp"] = df["overlap_end"] - df["overlap_start"]

    return df


def process_bed_chunk(bed_chunk, gtf, args):
    """
    Read in the chunked bed dataframe and the chromosome-specific GTF reference.
    Find the alignment/annotation overlaps using bioframe and assign a gene
    based on the results.

    :param args: object containing all supplied arguments
    :type args: class argparse.Namespace
    :return: dataframe of alignment intervals
    :rtype: pandas dataFrame
    """
    # The bed file has alignments and the chromosome has annotations,
    # so process the overlaps
    df_chunk = get_overlaps(bed_chunk, gtf)
    df_chunk["status"] = "Unknown"
    df_chunk = assign_status_low_mapq(df_chunk, args)
    df_chunk = assign_status_ambiguous_overlap(df_chunk)
    df_chunk = assign_status_no_features(df_chunk)
    df_chunk = find_largest_overlap(df_chunk)

    df_chunk = df_chunk[["read", "status", "score", "gene", "index_bed"]]

    df_chunk = df_chunk.reset_index(drop=True)
    df_chunk = df_chunk.drop(["index_bed"], axis=1)

    return df_chunk


def main(args):
    gtf = load_gtf(args)
    bed = load_bed(args)

    ext = args.output.split(".")[-1]
    if (bed.shape[0] > 0) & (gtf.shape[0] > 0):
        # Process alignment overlaps in chunks of <args.chunk_size> alignments
        n = int(np.ceil(bed.shape[0] / args.chunk_size))
        chunk_fns = []
        for i, bed_chunk in enumerate(np.array_split(bed, n)):
            df_chunk = process_bed_chunk(bed_chunk, gtf, args)

            fn = args.output.replace(ext, f"{i}.{ext}")
            chunk_fns.append(fn)
            df_chunk.to_csv(fn, sep="\t", index=False, header=False)

        # Concatenate chunked output files
        with open(args.output, "w") as f_out:
            for fn in chunk_fns:
                with open(fn) as f_in:
                    for line in f_in:
                        f_out.write(line)

        # Clean up chunked files
        [os.remove(fn) for fn in chunk_fns]

    else:
        # The bed file contained no alignments or the chromosome does not
        # have any annotations, so output empty file
        f = open(args.output, "w")
        f.close()


if __name__ == "__main__":
    args = parse_args()

    main(args)
