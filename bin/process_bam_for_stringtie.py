#!/usr/bin/env python
"""Preprocess bam ready for transcript assemble with stringtie.

Toggle the 16/0 bam bit to effectively invert the alignment.
This puts the alignemnt into forward orientation as tehy are in reverse
orientation after stranding.

Also split and write the new bams by barcode. Keep only the largest read for
each UMI.
TODO: Add a umi-guided consensus generation step.
"""
from collections import defaultdict
import sys

import pysam


def main():
    """Entry point."""
    bam_file = sys.argv[1]

    bc_alns = defaultdict(dict)

    with pysam.AlignmentFile(bam_file, "rb") as bam:

        for align in bam.fetch():
            bc = align.get_tag('CB')
            umi = align.get_tag('UB')
            aln_len = align.query_length
            flag = align.flag
            if flag == 16:
                align.flag = 0
            elif flag == 0:
                align.flag = 16
            else:
                raise ValueError('Only expecting bam flags of 16 or 0')
            # Keep only the UMI with the longest read
            try:
                bc_alns[bc][umi]
            except KeyError:
                bc_alns[bc][umi] = align
            else:
                if aln_len > bc_alns[bc][umi].query_length:
                    bc_alns[bc][umi] = align

    for bc_, umi_dicts in bc_alns.items():
        outfile = f'{bc_}.bam'
        with pysam.AlignmentFile(outfile, 'wb', template=bam) as bamout:
            for _, aln_ in umi_dicts.items():
                bamout.write(aln_)


if __name__ == "__main__":
    main()
