### Input Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| fastq | string | FASTQ files to use in the analysis. | This accepts one of three cases: (i) the path to a single FASTQ file; (ii) the path to a top-level directory containing FASTQ files; (iii) the path to a directory containing one level of sub-directories which in turn contain FASTQ files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| bam | string | BAM or unaligned BAM (uBAM) files to use in the analysis. | This accepts one of three cases: (i) the path to a single BAM file; (ii) the path to a top-level directory containing BAM files; (iii) the path to a directory containing one level of sub-directories which in turn contain BAM files. In the first and second case, a sample name can be supplied with `--sample`. In the last case, the data is assumed to be multiplexed with the names of the sub-directories as barcodes. In this case, a sample sheet can be provided with `--sample_sheet`. |  |
| ref_genome_dir | string | The path to the 10x reference directory | Human reference data can be downloaded from 10x [here](https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz). Instructions for preparing reference data can be found [here](https://www.10xgenomics.com/support/software/cell-ranger/tutorials/cr-tutorial-mr#overview) |  |
| kit | string | The 10x kit and version separated by a colon (eg: 3prime:v3) | 10x kits can be released with different versions, each requiring a specific whitelist that is looked-up by the workflow. If `single_cell_sample_sheet` is not defined, the 10x kit is applied to all samples. This parameter is ignored if `single_cell_sample_sheet` is supplied. |  |
| expected_cells | integer | Number of expected cells in the sample. | The number of expected cells. If `single_cell_sample_sheet` is not defined, `expected_cells` is applied to all samples. This parameter is ignored if `single_cell_sample_sheet` is supplied. |  |
| single_cell_sample_sheet | string | An optional CSV file used to assign library metadata per sample. If all samples have the same library metadata, this can be supplied instead by using the `--kit` and `--expected_cells` parameters. | Columns should be: [sample_id, kit, exp_cells]. This must not be confused with the MinKNOW sample_sheet. `sample_id` should correspond to `sample_name` which is defined either in the `sample_sheet`, given by the `sample` parameter (for single sample runs) or if no `sample_sheet` or `sample` is given, is derived from the folder name containing the FASTQ files. |  |
| full_length_only | boolean | Only process full length reads. | If set to true, only process reads or subreads that are classified as full length (read segments flanked by compatible adapters in the expected orientation). | True |
| min_read_qual | number | Specify read quality lower limit. | Any reads with a quality lower than this limit will not be included in the analysis. |  |


### Sample Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| sample_sheet | string | A CSV file used to map barcodes to sample aliases. The sample sheet can be provided when the input data is a directory containing sub-directories with FASTQ files. | The sample sheet is a CSV file with, minimally, columns named `barcode` and `alias`. Extra columns are allowed. A `type` column is required for certain workflows and should have the following values; `test_sample`, `positive_control`, `negative_control`, `no_template_control`. |  |
| sample | string | A single sample name for non-multiplexed data. Permissible if passing a single .fastq(.gz) file or directory of .fastq(.gz) files. |  |  |


### Output Options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| out_dir | string | Directory for output of all workflow results. |  | output |


### Advanced options

| Nextflow parameter name  | Type | Description | Help | Default |
|--------------------------|------|-------------|------|---------|
| kit_config | string | A file defining the configurations associated with the various supported 10x kits. | A CSV file is expected with the following headers [kit, barcode_length, umi_length]. If not specified, a default `kit_configs.csv` (found in the project directory root) will be used. This parameter does not typically need be changed. |  |
| threads | integer | Number of CPU threads to use in resource intensive processes. | The total CPU resource used by the workflow is constrained by the executor configuration. | 8 |
| fastq_chunk | integer | Sets the maximum number of reads per chunk for the initial processing of reads. | Controls batching of reads for processing. | 1000000 |
| barcode_adapter1_suff_length | integer | Suffix length of the read1 adapter to use in creating the probe sequence for identifying barcode/UMI bases. | For example, specifying 12 would mean that the last 12 bases of the specified read1 sequence will be included in the probe sequence. | 10 |
| barcode_min_quality | integer | Minimum allowed nucleotide-level quality score in the extracted/uncorrected barcode sequence. | Values equal or higher to this this will be considered 'high-quality' and used for generating the barcode whitelist. | 15 |
| barcode_max_ed | integer | Maximum allowable edit distance between uncorrected barcode and the best matching corrected barcode from the sample whitelist. | Barcodes are corrected by searching from a list of barcodes known to exist in the dataset. A maximum edit distance of 2 between query and whitelist barcode is recommended. | 2 |
| barcode_min_ed_diff | integer | Minimum allowable edit distance difference between whitelist barcode candidates. | If there is more than one candidate barcode found in the whitelist, the edit distance difference of the top hit and second best hits (in relation to the uncorrected barcode) must be at least this value to be able to assign a barcode. If the edit distance difference is less than this, it is assumed that barcode identity is amiguous, and the read is not tagged with a corrected barcode. | 2 |
| gene_assigns_minqv | integer | Minimum MAPQ score allowed for a read to be assigned to a gene. |  | 30 |
| matrix_min_genes | integer | Filter cells from the gene expression matrix if they contain fewer than <matrix_min_genes> genes. |  | 200 |
| matrix_min_cells | integer | Filter genes from the gene expression matrix that are observed in fewer than <matrix_min_cells> cells. |  | 3 |
| matrix_max_mito | integer | Filter cells from the gene expression matrix if more than <matrix_max_mito> percent of UMI counts come from mitochondrial genes. |  | 20 |
| matrix_norm_count | integer | Normalize expression matrix to <matrix_norm_count> counts per cell. |  | 10000 |
| genes_of_interest | string | File containing a list of gene symbols (one symbol per line) to annotate with expression values in the UMAP projections. If doing visium spatial analysis, these genes will be used to annotate the spatial plots.  |  |  |
| mito_prefix | string | Gene name prefix to identify for mitochondrial genes. | Parts of the workflow analyse mitochondrial genes separately. These genes are identified by searching for a gene name prefix. Human mitochondrial genes can be identified with prefix 'MT-' and mouse genes with prefix 'mt-'. If the reference genome contains data from multiple organisms with different nomenclature, multiple prefixes can be supplied like so: 'MT-,mt-' | MT- |
| umap_n_repeats | integer | Number of UMAP projection to repeat for each dataset. | The UMAP algorithm contains elements of randomness that can mislead users into seeing associations between cells that are not meaningful. It is recommended to view multiple plots generated with the same parameters and check that any observed structure is consistent across runs. | 3 |
| stringtie_opts | string | StringTie options for transcriptome assembly. | StringTie option string can be supplied at the command line as in this example: `--stringtie_opts="-c 5 -m 100 "`. StringTie options can be found here: http://ccb.jhu.edu/software/stringtie/index.shtml?t=manual. The default option (-c 2) ensures that only transcripts with a coverage of 2 or higher are included in the generated transcriptome | -c 2 |


