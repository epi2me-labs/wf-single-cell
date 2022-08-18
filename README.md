# Workflow single-cell

wf-single-cell is a research pipeline designed to identify the cell barcode
and UMI sequences present in nanopore sequencing reads generated from single-cell gene expression libraries. 

It was initially created as a Nextflow port of [Sockeye](https://github.com/nanoporetech/sockeye).




## Introduction

The following single-cell kits from 10x Genomics are currently supported:
- Chromium Single Cell [3ʹ gene expression](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html), versions 2 and 3
- Chromium Single Cell [5ʹ gene expression](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html), version 1
- Chromium Single Cell [Multiome (ATAC + GEX)](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html), version 1

Oxford Nanopore has developed a protocol for sequencing single-cell libraries from 10X, which can be found on the Nanopore Community [website](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/single-cell-transcriptomics-10x/v/sst_v9148_v111_revb_12jan2022).

The inputs to Sockeye are raw nanopore reads (FASTQ) generated from the sequencing
instrument and reference files that can be downloaded from [10X](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).
The pipeline outputs a gene x cell expression matrix, as well as a BAM file of
aligned reads tagged with cell barcode and UMI information.

Package dependencies
--------------------

The wf-single-cell pipeline makes use of the following dependencies.

- bedtools
- bioframe
- biopython
- editdistance
- matplotlib
- minimap2
- numpy
- pandas
- parasail-python
- pysam
- samtools
- scikit-learn
- seqkit
- tqdm
- umap-learn
- vsearch

Additionally, while no explicit dependency exists for the
[UMI-tools](https://github.com/CGATOxford/UMI-tools) package, the script
``bin/cluster_umis.py`` makes significant use of several functions from
the package. More detailed acknowledgements can be found in the source code.
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
nextflow run epi2me-labs/wf-template --help
```

to see the options for the workflow.

**Downloading reference data**
The pipeline requires access to reference data files that are packaged and freely available from 10x Genomics. For human samples, the GRCh38 packaged reference files can be downloaded using either curl or wget using:

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

or 
```
cd /PATH/TO/10X/DOWNLOADS
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

**Download demonstration data**
```
A small test dataset is provided for the purposes of testing the workflow software,
it consits of data from 10 cells and the 10x reference data for only chr22
it can be downloaded using:

wget -O test_data.tar.gz  \
  https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-single-cell/wf-single-cell-testdata-chr22.tar.gz \
  && tar -xzvf test_data.tar.gz
```

The single-cell-sample sheet defines the samples and 10x kits that were used to process them. 
An exampke sheet with one sample is:
```
run_id,kit_name,kit_version
run1,3prime,v3
```

The `run_id` must match a sample given in the `--fastq` option

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
## Useful links

* [nextflow](https://www.nextflow.io/)
* [docker](https://www.docker.com/products/docker-desktop)
* [conda](https://docs.conda.io/en/latest/miniconda.html)
