"""Test assign_barcodes."""
import pandas as pd
import pytest
from workflow_glue.create_matrix import cluster_dataframe


@pytest.fixture()
def umi_gene_df():
    """Make read tag and feature assignment DataFrames."""
    # Define 3 clusters of UMIs.
    # each entry contains (UMI, gene name and number of UMIs
    clusters = [

        # Cluster1 - 3 UMIS ###################
        # 'true' UMI
        ('AAAAAAAAAAAA', 'YFG1', 20),  # umi1
        # ED to umi1 = 2. n_true > (n_umi2 * 2) - 1
        ('ttAAAAAAAAAA', 'YFG1', 10),  # umi2

        # Cluster2 - single UMI ###############
        # ED to umi1 = 2, n_true < (n_umi2 * 2) - 1
        ('ggAAAAAAAAAA', 'YFG1', 15),  # umi3

        # Cluster3 - single UMI ###############
        # ED to umi1 = 3, n_true > (n_umi2 * 2) - 1
        ('AAAAAAAAAggg', 'YFG1', 10),  # umi4

        # Cluster4 - single UMI ###############
        # ED to umi1 = 1, n_true > (n_umi2 * 2) - 1,
        # but has a diffrent gene assignment to UMI 1
        ('cAAAAAAAAAAA', 'YFG2', 10),  # umi5

    ]

    # The actual dataframes used in the workflow will contain more columns,
    # but they are not used in the clustering process so are omitted for clarity.
    # CB is required, but can be any non '-' string
    header = ('read_id', 'UR', 'gene', 'CB')
    records = []
    read_num = 0

    for umi, gene, n_molecules in clusters:
        for _ in range(n_molecules):
            records.append((f'read_{read_num}', umi, gene, 'CB'))
            read_num += 1

    df = pd.DataFrame(
        records, columns=header).set_index(
        'read_id', drop=True)

    return df


def test_process_records(umi_gene_df):
    """Check that process_records is clustering and correcting UMIs appropriately."""
    cluster_dataframe(umi_gene_df, 1000)

    assert 'UB' in umi_gene_df
    # Check for the correct number of clusters
    assert umi_gene_df['UB'].nunique() == 4

    # Check that UMI2 is corrected to the 'true' UMI of cluster1
    assert all(
        umi_gene_df.loc[
            umi_gene_df['UR'] == 'ttAAAAAAAAAA'].loc[:, 'UB'] == 'AAAAAAAAAAAA')

    # Check that the rest of the UMIs map back to themselves
    # as they are all single-UMI clsuters
    df_no_clust1 = umi_gene_df[~umi_gene_df.UR.isin(['AAAAAAAAAAAA', 'ttAAAAAAAAAA'])]
    assert all(df_no_clust1.UR == df_no_clust1.UB)
