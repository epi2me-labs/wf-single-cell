
"""Make sparse snv x cell matrix.

TODO: We may want the expression matrix to use this code.
Or modify the expression matrix code to deal with genotype matrices
"""
import gzip
from pathlib import Path
import sys

import pandas as pd
import pysam


from .util import wf_parser  # noqa: ABS101


def argparser():
    """Create argument parser."""
    parser = wf_parser("MEX matrix")

    parser.add_argument(
        "vcf_in", help="VCF file")
    parser.add_argument(
        "mex_out_dir", help="MEX output directory")
    parser.add_argument(
        "--report_vars",
        help="List of variant ot report (chr_pos_ref_alt)", default=None)

    return parser


def main(args):
    """Write a matrix to disk in mtx format.

    :param args.matrix: Path to encoded genotype matrix file with encoding:
        hom ref: 0
        het: 1
        hom alt: 2
        no data: -1
    """
    # Full (non-sparse) matrix variants to write to a TSV file.
    for_report = []
    max_report_vars = 50
    report_vars_written = 0

    if args.report_vars:
        # Load the list of interesting variants
        with open(args.report_vars, 'r') as fh:
            report_vars = [line.strip() for line in fh]

    mex_folder = Path(args.mex_out_dir)
    mex_folder.mkdir(parents=True, exist_ok=True)

    vcf = pysam.VariantFile(args.vcf_in, threads=6)
    samples = list(vcf.header.samples)

    fhf = gzip.open(mex_folder / "features.tsv.gz", 'wt')

    n_rows = 0
    n_vars = 0

    with fhf as fh_feat:

        for i, record in enumerate(vcf.fetch()):
            n_rows += 1
            write_for_report = False
            if args.report_vars:
                if f"{record.chrom}_{record.pos}" in report_vars:
                    if report_vars_written < max_report_vars:
                        write_for_report = True
                        report_vars_written += 1
            else:
                # Just get the first n variants to show in the report
                if report_vars_written < max_report_vars:
                    write_for_report = True
                    report_vars_written += 1
            variant_id = f"{record.chrom}_{record.pos}_{record.ref}_{record.alts[0]}"
            fh_feat.write(variant_id + '\n')
            for j, sample in enumerate(samples):
                genotype = record.samples[sample]['GT']  # Get genotype
                # numerically-encode diploid genotype
                try:
                    gt_val = sum(allele for allele in genotype)
                except TypeError:
                    gt_val = -1  # Missing genotype
                if write_for_report:
                    for_report.append(
                        (variant_id, sample, gt_val))
                if gt_val == -1:
                    continue  # Skip missing genotypes from sparse matrix
                n_vars += 1
                # 1-based indexing for mtx format
                sys.stdout.write(f"{i + 1} {j + 1} {gt_val}\n")
    with gzip.open(mex_folder / "barcodes.tsv.gz", 'wt') as fh:
        for col in samples:
            fh.write(f"{col}\n")

    header = ((
            '%%MatrixMarket matrix coordinate integer general\n'
            '%metadata_json:'
            '{"software_version": "ont-single-cell","format_version": 2}\n'
            f'{n_rows} {len(samples)} {n_vars} \n')
        )
    with open('header.txt', 'w') as f:
        f.write(header)

    # Write full matrix of interesting variants to a TSV file
    df_top = pd.DataFrame.from_records(
        for_report, columns=['variant', 'sample', 'gt_val'])
    df_top = df_top.pivot(
        index='variant', columns='sample', values='gt_val').fillna(-1)
    df_top.to_csv("top_snvs.tsv", sep='\t', index=True)
