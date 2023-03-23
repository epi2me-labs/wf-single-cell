"""Test assign_barcodes."""
import tempfile

import pandas as pd
from workflow_glue.assign_features import (
    main
)


def test_main():
    """Test main."""
    # The read to stringtie transcrtipt dataframe from align_to_transcriptome
    df_query_transcript_rows = (
        ('read_1', 'ST001'),
        ('read_2', 'ST002'),
        ('read_3', 'ST003'),
        ('read_4', 'ST004')
    )
    df_stringtie = pd.DataFrame(df_query_transcript_rows)
    transcript_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_stringtie.to_csv(
        transcript_file, sep='\t', index=None, header=['read_id', 'qry_id'])

    # gffcompare tmap dataframe. Maps stringtie transcripts (qry_id)
    # to reference transcripts and reference gene IDs
    df_gffcompare_tmap_rows = (
        ('ST001', 'ref_tr_1', 'gene_id_1', '='),
        ('ST002', 'ref_tr_2', 'gene_id_2', '='),
        ('ST003', 'ref_tr_3', 'gene_id_3', '='),
        ('ST004', 'ref_tr_4', 'gene_id_4', 'i'),  # Contained solely withon an intron
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
        ('read_3', '40'),  # < mapq threxhold of 60
        ('read_4', '60')
    )
    df_tags = pd.DataFrame(
        df_tags_rows, columns=['read_id', 'mapq']
    )
    tags_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_tags.to_csv(tags_file, index=None, sep='\t')

    # GTF file, for converting gene_ids in the gffcomapre tmap file to gene_name./
    # Just a sibset of the gtf. We need transcript in pos [2] and gene_id and
    # transcript id are grepped
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
        query_transcript_read_assign = transcript_file
        gffcompare_tmap = gffcompare_file
        tags = tags_file
        gtf = gtf_file
        output = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False).name
        min_mapq = 60

    args = Args()
    main(args)

    result = pd.read_csv(args.output, sep='\t', index_col=0)

    assert result.at['read_1', 'gene'] == 'gene_name_1'
    assert result.at['read_1', 'transcript'] == 'ref_tr_1'

    assert result.at['read_3', 'gene'] == '-'
    assert result.at['read_3', 'transcript'] == '-'

    assert result.at['read_4', 'gene'] == 'gene_name_4'
    assert result.at['read_4', 'transcript'] == '-'
