"""Join tags from two TSV files based on read IDs."""
from contextlib import ExitStack
import sys
from .util import wf_parser  # noqa: ABS101


def argparser():
    """Parse the arguments."""
    parser = wf_parser(
        "Merge tag files.")
    parser.add_argument(
        "parent_file",
        help="Sorted tags TSV file containing all read IDs.")
    parser.add_argument(
        "child_file",
        help="Sorted tags TSV file containing a subset of read IDs from parent file.")
    parser.add_argument(
        "--output_file",
        help="Output file to write merged tags. If not provided, writes to stdout.")

    return parser


def main(args):
    """Join sorted parent and child files on read IDs.

    Child file may be a subset of parent.
    To avoid reading all into memory, read line by line. If there's a match, write the
    merged row.
    """
    with ExitStack() as stack:
        bufsize = 4*1024*1024  # 4 MiB buffer size
        fh_parent = stack.enter_context(
            open(args.parent_file, 'r', buffering=bufsize))
        fh_child = stack.enter_context(
            open(args.child_file, 'r', buffering=bufsize))
        if args.output_file:
            fh_output = stack.enter_context(
                open(args.output_file, 'w', buffering=bufsize))
        else:
            fh_output = sys.stdout
        parent_line = fh_parent.readline()
        child_line = fh_child.readline()

        while parent_line and child_line:
            # TODO: no need to split each time around
            parent_fields = parent_line.rstrip('\n').split('\t')
            child_fields = child_line.rstrip('\n').split('\t')

            parent_id = parent_fields[0]
            child_id = child_fields[0]

            if parent_id == child_id:
                # Match: join lines and write columns, don't duplicate the read ID
                merged_line = '\t'.join(parent_fields + child_fields[1:])
                fh_output.write(merged_line + '\n')
                parent_line = fh_parent.readline()
                child_line = fh_child.readline()

            elif parent_id < child_id:
                # parent is behind, so advance it to catch up.
                parent_line = fh_parent.readline()

            elif parent_id > child_id:
                # Child is behind. This can happen as some of the long reads may not
                # be in the tags file as we filter unmapped reads in the mapping process
                child_line = fh_child.readline()
