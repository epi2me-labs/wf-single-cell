"""Common code to be used across workflow scripts."""
import collections
import json

kit_adapters = {
    '3prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    'multiome': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    'visium': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'ATGTACTCTGCGTTGATACCACTGCTT'
    },
    '5prime': {
        'adapter1': 'CTACACGACGCTCTTCCGATCT',
        'adapter2': 'GTACTCTGCGTTGATACCACTGCTT'
    }
}

revcomp_map = str.maketrans("ACGTacgt", "TGCAtgca")


def rev_cmp(seq):
    """Reverse complement a DNA sequence."""
    return seq[::-1].translate(revcomp_map)


class StatsSummary(collections.Counter):
    """Summary dictionary for storing."""

    fields = {}  # subclasses should fill this in

    def __init__(self, *args, **kwargs):
        """Count some numbers."""
        self.update(*args, **kwargs)

    @classmethod
    def from_pandas(cls, df):
        """Create an instance from a pandas dataframe."""
        raise NotImplementedError("This method has not been implemented.")

    def to_dict(self):
        """Create dictionary with explicit zeroes."""
        return {k: self[k] for k in self}

    @classmethod
    def from_json(cls, fname):
        """Create and instance from a JSON file."""
        with open(fname, "r") as fh:
            data = json.load(fh)
        return cls(data)

    def to_json(self, fname):
        """Save to JSON."""
        with open(fname, "w") as fh:
            json.dump(self.to_dict(), fh, indent=4)
