"""Test assign_barcodes."""
import tempfile

import pandas as pd
from workflow_glue.assign_features import (
    main
)


def test_main():
    """Test main."""
    df_query_transcript_rows = (
        ('read_1', 'ST001'),
        ('read_2', 'ST002'),
        ('read_3', 'ST003'),
        ('read_4', 'ST004')
    )
    df_stringtie = pd.DataFrame(df_query_transcript_rows)
    transcript_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_stringtie.to_csv(transcript_file, sep='\t', index=None, header=None)

    df_gffcompare_tmap_rows = (
        ('ST001', 'YFT1', 'YFG1', '='),
        ('ST002', 'YFT2', 'YFG2', '='),
        ('ST003', 'YFT3', 'YFG3', '='),
        ('ST004', 'YFT4', 'YFG4', 'i'),  # Contained solely withon an intron
    )
    df_gffcompare_tmap = pd.DataFrame(
        df_gffcompare_tmap_rows, columns=[
            ['qry_id', 'ref_id', 'ref_gene_id', 'class_code']]
    )
    gffcompare_file = tempfile.NamedTemporaryFile('w', delete=False, suffix='.tsv').name
    df_gffcompare_tmap.to_csv(gffcompare_file, sep='\t', index=None)

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

    class Args:
        query_transcript_read_assign = transcript_file
        gffcompare_tmap = gffcompare_file
        tags = tags_file
        output = tempfile.NamedTemporaryFile('w', suffix='.tsv', delete=False).name
        min_mapq = 60

    args = Args()
    main(args)

    result = pd.read_csv(args.output, sep='\t', index_col=0)

    assert result.at['read_1', 'gene'] == 'YFG1'
    assert result.at['read_1', 'transcript'] == 'YFT1'

    assert result.at['read_3', 'gene'] == '-'
    assert result.at['read_3', 'transcript'] == '-'

    assert result.at['read_4', 'gene'] == 'YFG4'
    assert result.at['read_4', 'transcript'] == '-'
