## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[conda](https://docs.conda.io/en/latest/miniconda.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or conda is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-single-cell --help
```

to see the options for the workflow.


The main options are:
* `fastq`: A fastq file or directory containing fastq input files or directories of input files.
* `ref_genome_dir` The path to the 10x reference genome directory (see `Downloading reference data` below)
* `single_cell_sample_sheet`
This is used to associate the samples with 10x kits that were used to process them, and
is not to be confused with the optional `sample_sheet` from the basecaller. 
The run_id should correspond to sample_id which is defined either in the `sample_sheet`, given by the `sample`
parameter (for single sample runs) or if no `sample_sheet` or `sample` is given, is derived from the folder name containing the fastq files.

An example sheet with one sample is:
```
run_id,kit_name,kit_version
run1,3prime,v3
```

**Downloading reference data**
The pipeline requires access to reference data files that are packaged and freely available from 10x Genomics. For human samples, the GRCh38 packaged reference files can be downloaded using either curl or wget using:

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

or 
```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

**Download demonstration data**
```
A small dataset is provided for the purposes of testing the workflow
It consits of data from just 10 cells and the 10x reference data for only chr22.
It can be downloaded using:

wget -O test_data.tar.gz  \
  https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-single-cell/wf-single-cell-testdata-chr22.tar.gz \
  && tar -xzvf test_data.tar.gz
```

The workflow can be run with the demonstration data using:

```
OUTPUT=output
nextflow run epi2me-labs/wf-single-cell \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq wf-single-cell-testdata-chr22/fastq \
    --single_cell_sample_sheet wf-single-cell-testdata-chr22/samples.csv \
    --ref_genome_dir wf-single-cell-testdata-chr22/refdata-gex-GRCh38-2020-A-chr22 \
    --out_dir ${OUTPUT}
```

The output of the pipeline will be found in `./output` for the above
example. 

**Workflow outputs**

The pipeline output will be written to a directory defined by ``--out_dir``. 
For each sample specifed in the `single_cell-sample_sheet`  an output folder is created containing the results:


The most useful outputs of the pipeline are likely:

* ``configs.stats.json``: provides a summary of sequencing statistics and observed read configurations, such as

  - ``n_reads``: number of total reads in the input fastq(s)
  - ``rl_mean``: mean read length
  - ``n_fl``: total number of reads with the read1-->TSO or TSO'-->read1' adapter configuration (i.e. full-length reads)
  - ``n_plus``: number of reads with the read1-->TSO configuration
  - ``n_minus``: number of reads with the TSO'-->read1' configuration

* ``tagged.sorted.bam``: BAM file of alignments to the reference where each alignment contains the following sequence tags

  - CB: corrected cell barcode sequence
  - CR: uncorrected cell barcode sequence
  - CY: Phred quality scores of the uncorrected cell barcode sequence
  - UB: corrected UMI sequence
  - UR: uncorrected UMI sequence
  - UY: Phred quality scores of the uncorrected UMI sequence

* ``gene_expression.processed.tsv``: TSV containing the gene (rows) x cell (columns) expression matrix, processed and normalized according to the folloing parameters:

  - ``matrix_min_genes``: cells with fewer than this number of expressed genes will be removed
  - ``matrix_min_cells``: genes present in fewer than this number of cells will be removed
  - ``matrix_max_mito``: cells with more than this percentage of counts belonging to mitochondrial genes will be removed
  - ``matrix_norm_count``: normalize all cells to this number of total counts per cell
