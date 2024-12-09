"""Add tags from a text file to a BAM."""
import csv
from dataclasses import dataclass
import itertools
from pathlib import Path

import pysam
from .util import get_named_logger, wf_parser  # noqa: ABS101

logger = get_named_logger("TagBAMs")


BAM_TAGS = {
    "corrected_barcode": "CB",
    "uncorrected_barcode": "CR",
    "quality_barcode": "CY",
    "corrected_umi": "UB",
    "uncorrected_umi": "UR",
    "quality_umi": "UY",
    "gene": "GN",
    "transcript": "TR"
}


def argparser():
    """Create argument parser."""
    parser = wf_parser("tag_bams")

    parser.add_argument(
        "in_bam", type=Path,
        help="BAM file for tagging")

    parser.add_argument(
        "out_bam", type=Path,
        help="Path for tagged output BAM")

    parser.add_argument(
        "tags", type=Path,
        help="Read tags TSV")

    parser.add_argument(
        "--threads", default=2, type=int,
        help="Number of threads used for BAM reading/writing.")
    return parser


# The use of a dataclass here is primarily to reduce memory:
#  chr1, 9517964 reads. dict: 11.6 GB, class: 8.5 GB
# The overhead in creating instances of the class is small
# compared to the BAM writing time, and access is similarly
# fast enough.

@dataclass
class Tags:
    """Storing tag data for a read."""

    CB = None
    CR = None
    CY = None
    UB = None
    UR = None
    UY = None
    GN = None
    TR = None
    chrom = None

    @classmethod
    def from_dict(cls, d):
        """Create instance from a dictionary."""
        self = cls()
        for k in BAM_TAGS.values():
            setattr(self, k, d[k])
        setattr(self, "chrom", d["chr"])
        return self


class TagStore:
    """Proxy to tag files for retrieving per-read tag information."""

    def __init__(self, tags, bam=None):
        """Initialize an instance."""
        self._tags = None
        self._cur = None
        self._single = False
        if tags.is_file():
            self._single = True
            self._tags = self._read_file(tags)
        elif tags.is_dir():
            if bam is None:
                raise ValueError("`bam` should be provided when `tags` is a directory.")
            tags = tags.glob("*.tsv")
            self._index = dict()
            for fname in tags:
                d = self._read_file(fname, nrows=10)
                try:
                    chrom = getattr(next(iter(d.values())), "chrom")
                    self._index[chrom] = fname
                except StopIteration:
                    logger.warn(f"{fname} appears empty.")
                else:
                    logger.info(f"{fname} contains tags for reference: {chrom}.")
        else:
            raise ValueError(
                "`tags` should be a tags file or directory containing such files.")

    def populate(self, rname):
        """Populate the proxy for a given reference."""
        if not self._single:
            self._cur = rname
            self._tags = self._read_file(self._index[self._cur])

    def _read_file(self, fname, nrows=None, cols=None):
        """Read a tags file."""
        # note: this is actually around 50% faster than:
        #       pd.read_csv().to_dict(orient="index")
        # first find and rename the fields to be tags rather than human names
        fields = None
        with open(fname) as csvfile:
            iterator = csv.DictReader(csvfile, delimiter="\t")
            fields = iterator.fieldnames
            for k, v in BAM_TAGS.items():
                for i in range(len(fields)):
                    if fields[i] == k:
                        fields[i] = v
        # now parse the file
        with open(fname) as csvfile:
            iterator = csv.DictReader(csvfile, delimiter="\t", fieldnames=fields)
            next(iterator)  # setting fieldnames doesn't read header
            if nrows is not None:
                iterator = itertools.islice(iterator, nrows)
            data = {d["read_id"]: Tags.from_dict(d) for d in iterator}
        return data

    def __getitem__(self, read_data):
        """Retrieve tags for a read."""
        read_id, chrom = read_data
        try:
            data = self._tags[read_id]
        except KeyError:
            exp = KeyError(f"Read '{read_id}' not found in tag data.")
            if chrom != self._cur:
                self._cur = chrom
                self._tags = self._read_file(self._index[self._cur])
                try:
                    data = self._tags[read_id]
                except KeyError:
                    raise exp
            else:
                raise exp
        return data


def add_tags(tags, in_bam, out_bam, threads):
    """Add all the required tags to the BAM file."""
    store = TagStore(tags, bam=in_bam)

    skipped = 0
    written = 0
    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in:
        with pysam.AlignmentFile(
                out_bam, "wb", template=bam_in, threads=threads) as bam_out:
            for ref in bam_in.references:
                logger.info(f"Processing reads from reference: {ref}.")
                try:
                    store.populate(ref)
                except KeyError:
                    logger.warn(f"Could not find tag data for reference: {ref}.")
                    continue
                logger.info("Tagging reads.")
                for align in bam_in.fetch(ref):
                    read_id = align.query_name
                    try:
                        row = store._tags[read_id]
                    except KeyError:
                        skipped += 1
                        continue  # don't write reads without tags
                    else:
                        written += 1
                        for tag in BAM_TAGS.values():
                            align.set_tag(tag, getattr(row, tag), value_type="Z")
                        bam_out.write(align)
    total = skipped + written
    written_pct = 0
    skipped_pct = 0
    if total > 0:
        written_pct = 100 * written / total
        skipped_pct = 100 * skipped / total
    logger.info(
        f"Written: {written} ({written_pct:0.2f}%). "
        f"Skipped: {skipped} ({skipped_pct:0.2f}%).")


def main(args):
    """Entry point."""
    add_tags(args.tags, args.in_bam, args.out_bam, args.threads)
