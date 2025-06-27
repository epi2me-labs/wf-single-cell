"""Test the join_tags script."""
from unittest.mock import Mock

import pandas as pd
import pytest
from workflow_glue.join_tags import main


@pytest.mark.parametrize(
    "child_content, expected_reads",
    [
        # Valid subset — should pass
        (
            "read_id\tchild_umi\n"
            "read_1\tumi1\n"
            "read_3\tumi3\n"
            "read_5\tumi5\n"
            "read_6\tumi6\n",
            ['read_1', 'read_3', 'read_5', 'read_6']
        )
    ]
)
def test_main_parametrized(tmp_path, child_content, expected_reads):
    """Test the main function."""
    parent = tmp_path / 'parent_input.tsv'
    child = tmp_path / 'child_input.tsv'
    output = tmp_path / 'output.tsv'

    # Write the parent input (shared by both cases)
    parent.write_text(
        "read_id\tparent_umi\n"
        "read_1\tumi1\n"
        "read_2\tumi2\n"
        "read_3\tumi3\n"
        "read_4\tumi4\n"
        "read_5\tumi5\n"
        "read_6\tumi6\n"
        "read_7\tumi7\n"
    )

    # Write the child input
    child.write_text(child_content)

    args = Mock()
    args.parent_file = parent
    args.child_file = child
    args.output = output

    main(args)
    merged = pd.read_csv(output, sep='\t')
    assert list(merged['read_id']) == expected_reads
