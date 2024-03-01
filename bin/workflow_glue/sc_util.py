"""Common code to be used across workflow scripts."""

kit_adapters = {
    '3prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    'multiome': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    '5prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'GTACTCTGCGTTGATACCACTGCTT'
    }
}


def rev_cmp(seq):
    """Reverse complement a DNA sequence."""
    revcomp_map = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq[::-1].translate(revcomp_map)
