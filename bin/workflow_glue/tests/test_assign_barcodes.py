"""Test assign_barcodes."""
import tempfile

import pandas as pd
from pytest import fixture
from workflow_glue.assign_barcodes import (
    calc_ed_with_whitelist, process_records
)


@fixture
def whitelist():
    """Make a small barcode whitelist."""
    return (
        'AAAAAAAAAAAAAAAA',
        'ttAAAAAAAAAAAAAA',
        'ttttAAAAAAAAAAAA',
        'AAAAAAAAAAAAAccc',
        'AAAggggAAAAAAAAA')


def test_calc_ed_with_whitelist(whitelist):
    """Test edit distance calculation."""
    bc1 = 'AAAAAAAAAAAAAAAA'
    bc_match, bc_match_ed, next_match_diff = \
        calc_ed_with_whitelist(bc1, whitelist)
    # This BC should match the first BC in the whitelist with an ED of 0
    assert bc_match == 'AAAAAAAAAAAAAAAA'
    assert bc_match_ed == 0

    # AAAAAAAAAAAAAAAA -> 0
    # ttAAAAAAAAAAAAAA -> 2
    # Check ED difference between top match and second-top match
    assert next_match_diff == 2

    # An uncorrected BC with a nearest match ED of 7 (cutoff = 6)
    # will return a string of 'X's and an ED the size of the BC.
    bc2 = 'AAAAAAAAAggggggg'
    bc_match, bc_match_ed, next_match_diff = \
        calc_ed_with_whitelist(bc2, whitelist)
    assert bc_match == 'X' * 16
    assert bc_match_ed == 16
    assert next_match_diff == 16


def test_process_records(whitelist):
    """Test process_records.

    Check if barcodes are corrected and enumerted appropriately.
    """
    # Build some uncorrectred barcodes.
    # The minimum required columns are read_id and CR (uncorrected barcode)
    header = ('read_id', 'CR')
    rows = [
        # These should be corrected to AAAAAAAAAAAAAAAA
        ('read1', 'AAAAAAAAAAAAAAAA'),
        ('read2', 'AAAAcAAAcAAAAAAA'),
        # These will not be corrected as the edit distance difference between
        # the top match and the second-top match is < 2
        ('read3', 'tAAAAAAAAAAAAAAA'),
        ('read4', 'AtAAAAAAAAAAAAAA'),
    ]

    tags = pd.DataFrame(rows, columns=header).set_index('read_id', drop=True)
    tags_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv')
    tags.to_csv(tags_file.name, sep='\t')

    barcode_length = 16
    max_ed = 16
    min_ed_diff = 2
    barcode_counter, result_tags_df = \
        process_records(tags_file.name, whitelist, barcode_length, max_ed, min_ed_diff)

    # Just the single corrected barcode should be present: AAAAAAAAAAAAAAAA
    assert len(barcode_counter) == 1
    assert barcode_counter['AAAAAAAAAAAAAAAA'] == 2

    assert result_tags_df.loc['read1', 'CB'] == 'AAAAAAAAAAAAAAAA'
    assert result_tags_df.loc['read2', 'CB'] == 'AAAAAAAAAAAAAAAA'
    assert result_tags_df.loc['read3', 'CB'] == '-'
    assert result_tags_df.loc['read4', 'CB'] == '-'
