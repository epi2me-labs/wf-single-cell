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
        help="Path for the tagged output BAM",
        default='-')

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
                    logger.warning(f"{fname} appears empty.")
                else:
                    logger.info(f"{fname} contains tags for reference: {chrom}.")
        else:
            raise ValueError(
                "`tags` should be a tags file or directory containing such files.")

    def populate(self, rname):
        """Populate the proxy for a given reference."""
        if not self._single:
            self._cur = rname
            try:
                self._tags = self._read_file(self._index[self._cur])
            except KeyError:
                self._tags = {}

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

    n_skipped = 0
    n_prim_written = 0
    n_supp_written = 0

    threads = max(1, threads // 2)

    b1 = pysam.AlignmentFile(in_bam, "rb", threads=threads)
    # Uncompressed BAM if writing to stdout
    write_mode = 'wbu' if out_bam == "-" else 'wb'
    b2 = pysam.AlignmentFile(out_bam, write_mode, template=b1, threads=threads)

    sup_tags = {}

    with b1 as bam_in, b2 as bam_out:

        for ref in bam_in.references:
            logger.info(f"Processing reads from reference: {ref}.")
            store.populate(ref)

            logger.info("Tagging reads.")
            for align in bam_in.fetch(ref):
                read_id = align.query_name

                if align.is_supplementary:
                    # No tagging of supp. reads yet, that's done later
                    continue
                try:
                    row = store._tags[read_id]
                except KeyError:
                    # No tags for this read
                    n_skipped += 1
                    continue

                # If primary has tags and a supp record, save the ID and tags
                # for later tagging.
                try:
                    if align.get_tag("SA"):
                        tags = []
                        for tag in BAM_TAGS.values():
                            tags.append(getattr(row, tag))
                        sup_tags[read_id] = tags
                except KeyError:
                    pass  # No SA tag
                n_prim_written += 1
                for tag in BAM_TAGS.values():
                    align.set_tag(tag, getattr(row, tag), value_type="Z")
                bam_out.write(align)

        bam_in.reset()  # reset the iterator, this time tagging supp. records
        for record in bam_in.fetch():
            if not record.is_supplementary:
                continue
            read_id = record.query_name
            try:
                tags = sup_tags[read_id]
            except KeyError:
                # No tags for this read
                continue
            for tag, val in zip(BAM_TAGS.values(), tags):
                record.set_tag(tag, val, value_type="Z")
            n_supp_written += 1
            bam_out.write(record)

        total = n_skipped + n_prim_written
        written_pct = 0
        skipped_pct = 0
        if total > 0:
            written_pct = 100 * n_prim_written / total
            skipped_pct = 100 * n_skipped / total
        logger.info(
            f"Tagged primary records: {n_prim_written} ({written_pct:0.2f}%).\n"
            f"Skipped primary records: {n_skipped} ({skipped_pct:0.2f}%).\n"
            f"Tagged suppl. records: {n_supp_written}.")


def main(args):
    """Entry point."""
    add_tags(
        args.tags, args.in_bam, args.out_bam, args.threads)
