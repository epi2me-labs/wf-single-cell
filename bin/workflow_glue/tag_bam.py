"""Tag BAM with workflow-derived information.

Tags files are TSV files containing read_id to tags mappings (such as barcodes, UMIs,
assigned features). We iterate over the BAM file by chromosome, loading the tags for
each chromosome individually to avoid holding all tags in memory
at once. The records are tagged and output to a tagged BAM. This process only tags
primary records, or supplementary records that are on the same chromosome as their
primary record.

To tag supplementary records, there is another tags file input that contains
read_id to tag mappings for all supplementary records, formatted identically to the
primary tags file.
During tagging, these are all loaded into memory regardless of chromosome they map to.
This allows supplementary records that map to a different chromosome than their
primary alignment to be properly tagged.
"""
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
        "sa_tags", type=Path,
        help="Read supplementary tags TSV")

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

    def __init__(self, tags, bam=None, sa_tags=None):
        """Initialize an instance."""
        self._sa_tags = self._load_supplementary_tags(sa_tags)
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
                    logger.warning(f"{fname} appears empty.")
                else:
                    logger.info(f"{fname} contains tags for reference: {chrom}.")
        else:
            raise ValueError(
                "`tags` should be a tags file or directory containing such files.")

    def _load_supplementary_tags(self, sa_tags):
        # Load tag info for all reads with one or more suppl. records.
        # These are added to self._tags later regardless of chr mapping.
        # This enures that supplementary records are tagged even if on a different
        # chr to primary record.
        sa_tags_files = []
        if sa_tags is not None:
            sa_tags_files = sa_tags.glob("*.tsv")
        sa_tags = {}
        for fname in sa_tags_files:
            sa_tags.update(self._read_file(fname))
        return sa_tags

    def populate(self, rname):
        """Populate the proxy for a given reference."""
        if not self._single:
            self._cur = rname
            try:
                self._tags = self._read_file(self._index[self._cur])
            except KeyError:
                # No primary records for this chr, but there may be suppl records
                self._tags = {}
            self._tags.update(self._sa_tags)

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


def add_tags(tags, sa_tags, in_bam, out_bam, threads):
    """Add all the required tags to the BAM file."""
    store = TagStore(tags, bam=in_bam, sa_tags=sa_tags)

    skipped = 0
    written = 0
    with pysam.AlignmentFile(in_bam, "rb", threads=threads) as bam_in:
        with pysam.AlignmentFile(
                out_bam, "wb", template=bam_in, threads=threads) as bam_out:
            for ref in bam_in.references:
                logger.info(f"Processing reads from reference: {ref}.")
                # There may be no primary records for this reference, but we'll process
                # it in case there are any supplementary records
                store.populate(ref)
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
    add_tags(args.tags, args.sa_tags, args.in_bam, args.out_bam, args.threads)
