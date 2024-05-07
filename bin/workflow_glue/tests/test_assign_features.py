"""Test assign_barcodes."""
import subprocess
import tempfile

import pandas as pd
from workflow_glue.assign_features import (
    main
)


def getbam():
    """Create a synthetic bam file containing adapter and barcode sequences."""
    # two refrence transcripts in the BAM
    # One is 70 and one is 100bp long
    header = (
        "@SQ	SN:ST001	LN:2000\n"
        "@SQ	SN:ST002	LN:1000\n")

    # Define the transcripts. Normally we wouldn't have transcripts this samll as
    # they would not pass the size threadhold in the stringtie step, but they work for
    # out purposes here

    alns = [
        # Uniquely-mapped read
        ['read_1', 200, 100, 500, 60, 'ST001', '500M'],
        # read_2 mapped to two locations. ST002 should be assigned to read 2
        # as it's higher AS and read and transcript cov > 0.4
        ['read_2', 150, 100, 400, 1, 'ST001', '100H400M'],
        ['read_2', 200, 500, 500, 60, 'ST002', '500M'],
        # read_3 maps to two locations. The second alignemnt has higher AS score
        # but will not be assigned as reference coverage is < 0.4
        ['read_3', 150, 500, 100, 1, 'ST001', '50M400H50M'],
        ['read_3', 200, 500, 150, 60, 'ST001', '200H150M200H']
    ]

    sam = header
    for align in alns:
        qname, a_score, start, seqlen, mapq, rname, cigar = align
        # Make a sam file containing the read and a quality qscore of 60.
        sam += (
            f"{qname}\t0\t{rname}\t{start}\t{mapq}\t{cigar}\t*\t0\t0\t"
            f"{'A' * seqlen}\t{'?' * seqlen}\tAS:i:{a_score}\n"
        )

    # Write out a test BAM
    with tempfile.NamedTemporaryFile(
            mode='w', suffix='.sam', delete=False) as fh_sam:
        fh_sam.write(sam)
        sam_file = fh_sam.name

    bam = 'align.bam'
    subprocess.check_output(['samtools', 'view', sam_file, '-o', 'align.bam'])

    return bam


def test_main():
    """Test main."""
    # gffcompare tmap dataframe. Maps stringtie transcripts (qry_id)
    # to reference transcripts and reference gene IDs
    df_gffcompare_tmap_rows = (
        ('ST001', 'ref_tr_1', 'gene_id_1', '='),
        ('ST002', 'ref_tr_2', 'gene_id_2', '='),
        ('ST003', 'ref_tr_3', 'gene_id_3', '=')
    )
    df_gffcompare_tmap = pd.DataFrame(
        df_gffcompare_tmap_rows, columns=[
            ['qry_id', 'ref_id', 'ref_gene_id', 'class_code']]
    )
    gffcompare_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_gffcompare_tmap.to_csv(gffcompare_file, sep='\t', index=None)

    # All we want from tags is the mapq alignment score
    df_tags_rows = (
        ('read_1', '60'),
        ('read_2', '60'),
        ('read_3', '30'),
    )
    df_tags = pd.DataFrame(
        df_tags_rows, columns=['read_id', 'mapq']
    )
    tags_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_tags.to_csv(tags_file, index=None, sep='\t')

    # GTF file is used for mapping transcript id (from the gffcompare tmap file) to
    # gene name.
    # Here is just a subset of the gtf. We need 'transcript' in pos [2]. gene_name and
    # transcript_id are grepped
    gtf_str = (
        'chr1\tHAVANA\ttranscript\tgene_name "gene_name_1";transcript_id "ref_tr_1";\n',
        'chr1\tHAVANA\ttranscript\tgene_name "gene_name_2";transcript_id "ref_tr_2";\n',
        'chr1\tHAVANA\ttranscript\tgene_name "gene_name_3";transcript_id "ref_tr_3";\n',
        'chr1\tHAVANA\ttranscript\tgene_name "gene_name_4";transcript_id "ref_tr_4";',
    )

    with tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv') as fh:
        fh.writelines(gtf_str)
        gtf_file = fh.name

    class Args:
        transcriptome_bam = getbam()
        gffcompare_tmap = gffcompare_file
        tags = tags_file
        gtf = gtf_file
        output = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False).name
        min_mapq = 30
        min_tr_coverage = 0.4
        min_read_coverage = 0.4
        chunksize = 1

    args = Args()
    main(args)

    result = pd.read_csv(args.output, sep='\t', index_col=0)

    assert result.at['read_1', 'gene'] == 'gene_name_1'
    assert result.at['read_1', 'transcript'] == 'ref_tr_1'

    assert result.at['read_2', 'gene'] == 'gene_name_2'
    assert result.at['read_2', 'transcript'] == 'ref_tr_2'

    assert result.at['read_3', 'gene'] == 'gene_name_1'
    assert result.at['read_3', 'transcript'] == '-'
