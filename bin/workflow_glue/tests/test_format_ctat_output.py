"""Tests for the format_ctat_output module."""

from workflow_glue.format_ctat_output import load_fusion_data


def test_load_fusion_data_valid_file(tmp_path):
    """Test load_fusion_data with a valid fusion file."""
    # Create a temporary valid fusion file
    fusion_file = tmp_path / "fusions.tsv"
    fusion_file.write_text(
        "#FusionName\tLR_accessions\tLeftGene\tRightGene\tLeftBreakpoint"
        "\tRightBreakpoint\tSpliceType\n"

        "Fusion1\tread1,read2\tGeneA\tGeneB\tchr1:100\tchr2:200\tSpliceA\n"
        "Fusion2\tread3\tGeneC\tGeneD\tchr3:300\tchr4:400\tSpliceB\n"
    )

    fusion_dict = load_fusion_data(fusion_file)

    assert fusion_dict is not None
    assert len(fusion_dict) == 3  # 3 unique reads


def test_load_fusion_data_empty_file(tmp_path):
    """Test load_fusion_data with an empty fusion file."""
    fusion_file = tmp_path / "empty.tsv"
    fusion_file.write_text("")

    fusion_dict = load_fusion_data(fusion_file)

    assert fusion_dict is None


def test_load_fusion_data_no_entries(tmp_path):
    """Test load_fusion_data with a fusion file containing no valid entries."""
    # Create a fusion file with no valid entries
    fusion_file = tmp_path / "no_entries_fusion.tsv"
    fusion_file.write_text(
        "#FusionName\tLR_accessions\tLeftGene\tRightGene\tLeftBreakpoint"
        "\tRightBreakpoint\tSpliceType\n"
    )

    fusion_dict = load_fusion_data(fusion_file)

    assert fusion_dict is None


def test_load_fusion_data_duplicate_reads(tmp_path):
    """Test load_fusion_data with duplicate read IDs."""
    # Create a fusion file with duplicate read IDs
    fusion_file = tmp_path / "duplicate_reads_fusion.tsv"
    fusion_file.write_text(
        "#FusionName\tLR_accessions\tLeftGene\tRightGene\tLeftBreakpoint"
        "\tRightBreakpoint\tSpliceType\n"

        "Fusion1\tread1,read1\tGeneA\tGeneB\tchr1:100\tchr2:200\tSpliceA\n"
    )

    fusion_dict = load_fusion_data(fusion_file)

    assert fusion_dict is not None
    assert len(fusion_dict) == 1  # Only 1 unique read
    assert len(fusion_dict["read1"]) == 2  # 2 entries for the same read ID
