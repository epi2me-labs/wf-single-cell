"""Test assign_barcodes."""
from collections import Counter
import tempfile

import pandas as pd
import pytest
from workflow_glue.assign_barcodes import (
    determine_barcode, process_records
)


@pytest.fixture
def allowed_barcodes():
    """Make a small barcode whitelist."""
    return set(
        ('AAAAAAAAAAAAAAAA',
         'ttAAAAAAAAAAAAAA',
         'ttttAAAAAAAAAAAA',
         'AAAAAAAAAAAAAccc',
         'AAAggggAAAAAAAAA'))


def test_calc_ed_with_allowed_barcodes(allowed_barcodes):
    """Test edit distance calculation."""
    bc1 = 'AAAAAAAAAAAAAAAA'
    bc_match = determine_barcode(
        bc1, list(allowed_barcodes), allowed_barcodes, 2, 2, Counter())
    assert bc_match == 'AAAAAAAAAAAAAAAA'

    # An uncorrected BC with a nearest match ED of 7 (cutoff = 6)
    # return no match
    bc2 = 'AAAAAAAAAggggggg'
    bc_match = determine_barcode(
        bc2, list(allowed_barcodes), allowed_barcodes, 2, 2, Counter())
    assert bc_match == '-'


@pytest.mark.parametrize("use_kmer_index", [False, True])
def test_process_records(allowed_barcodes, use_kmer_index):
    """Test process_records.

    Check if barcodes are corrected and enumerted appropriately.
    """
    # Build some uncorrectred barcodes.
    # The columns used in this test are read_id and CR (uncorrected barcode). The other
    # Columns can be any value for now
    header = ('read_id', 'CR', 'CY', 'UR', 'UY', 'chr', 'start', 'end', 'mapq')
    rows = [
        # This should 100% match to whitelist
        ('read1', 'AAAAAAAAAAAAAAAA', 'qual', 'umi', 'qual' 'chr', 0, 100, 20),
        # This should be corrected to AAAAAAAAAAAAAAAA
        ('read2', 'AAAAcAAAcAAAAAAA', 'qual', 'umi', 'qual', 'chr', 0, 100, 20),
        # These will not be corrected as the edit distance difference between
        # the top match and the second-top match is < 2
        ('read3', 'tAAAAAAAAAAAAAAA', 'qual', 'umi', 'qual', 'chr', 0, 100, 20),
        ('read4', 'AtAAAAAAAAAAAAAA', 'qual', 'umi', 'qual', 'chr', 0, 100, 20),
        # No matches to whitelist
        ('read5', 'GGGGGGGGGGGGGGGG', 'qual', 'umi', 'qual', 'chr', 0, 100, 20),
        ('read6', 'GCGCGCGCGCGCGCGC', 'qual', 'umi', 'qual', 'chr', 0, 100, 20),
    ]

    tags = pd.DataFrame(rows, columns=header).set_index('read_id', drop=True)
    tags_file = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv')
    tags.to_csv(tags_file.name, sep='\t')
    tags_output = tempfile.NamedTemporaryFile(mode='w', suffix='.tsv')

    max_ed = 2
    min_ed_diff = 2
    barcode_counter, reasons_counter = process_records(
        tags_file.name, allowed_barcodes,
        max_ed, min_ed_diff, tags_output.name,
        use_kmer_index=use_kmer_index)

    result_tags_df = pd.read_csv(tags_output.name, sep='\t', index_col=0)

    # Just the single corrected barcode should be present: AAAAAAAAAAAAAAAA
    assert len(barcode_counter) == 1
    assert barcode_counter['AAAAAAAAAAAAAAAA'] == 2

    assert result_tags_df.loc['read1', 'CB'] == 'AAAAAAAAAAAAAAAA'
    assert result_tags_df.loc['read2', 'CB'] == 'AAAAAAAAAAAAAAAA'
    assert result_tags_df.loc['read3', 'CB'] == '-'
    assert result_tags_df.loc['read4', 'CB'] == '-'
    assert result_tags_df.loc['read5', 'CB'] == '-'
    assert result_tags_df.loc['read6', 'CB'] == '-'

    assert dict(reasons_counter) == \
           {
            'bc_shortlist_exact_match': 1,
            'bc_corrected': 1,
            'bc_no_shortlist_match': 2,
            'bc_shortlist_multiple_hits': 2
           }
