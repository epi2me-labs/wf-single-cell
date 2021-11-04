#!/usr/bin/env python
"""Script to check that sample sheet is well-formatted."""
import argparse
import sys

import pandas as pd

"""
Example minknow sample sheet
https://informatics-spec.git.oxfordnanolabs.local/minknow-specification/develop/sample_sheets.html
flow_cell_id,kit,sample_id,experiment_id,barcode,alias,type
FA026858,SQK-RBK004,barcoding_run,sequencing_20200522,barcode01,patient_id_5,test_sample
FA026858,SQK-RBK004,barcoding_run,sequencing_20200522,barcode02,patient_id_6,test_sample
FA026858,SQK-RBK004,barcoding_run,sequencing_20200522,barcode03,patient_id_7,test_sample
"""


def main():
    """Run entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument('sample_sheet')
    parser.add_argument('output')
    args = parser.parse_args()

    try:
        samples = pd.read_csv(args.sample_sheet, sep=None)
        if 'alias' in samples.columns:
            if 'sample_id' in samples.columns:
                sys.stderr.write(
                    "Warning: sample sheet contains both 'alias' and "
                    'sample_id, using the former.')
            samples['sample_id'] = samples['alias']
        if 'barcode' not in samples.columns \
                or 'sample_id' not in samples.columns:
            raise IOError()
    except Exception:
        raise IOError(
            "Could not parse sample sheet, it must contain two columns "
            "named 'barcode' and 'sample_id' or 'alias'.")
    # check duplicates
    dup_bc = samples['barcode'].duplicated()
    dup_sample = samples['sample_id'].duplicated()
    if any(dup_bc) or any(dup_sample):
        raise IOError(
            "Sample sheet contains duplicate values.")
    samples.to_csv(args.output, sep=",", index=False)


if __name__ == '__main__':
    main()
