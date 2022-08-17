## Introduction

The following single-cell kits from 10x Genomics are currently supported:
- Chromium Single Cell `3ʹ gene expression <https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html>`_, versions 2 and 3
- Chromium Single Cell `5ʹ gene expression <https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html>`_, version 1
- Chromium Single Cell `Multiome (ATAC + GEX) <https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html>`_, version 1

Oxford Nanopore has developed a protocol for sequencing single-cell libraries from 10X, which can be found on the Nanopore Community `website <https://community.nanoporetech.com/docs/prepare/library_prep_protocols/single-cell-transcriptomics-10x/v/sst_v9148_v111_revb_12jan2022>`_.

The inputs to Sockeye are raw nanopore reads (FASTQ) generated from the sequencing
instrument and reference files that can be downloaded from `10X
<https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest>`_.
The pipeline outputs a gene x cell expression matrix, as well as a BAM file of
aligned reads tagged with cell barcode and UMI information.

Prerequisites
-------------

``conda`` must be installed in order to create the base environment where the
Sockeye snakemake pipeline will run. Installation instructions can be found in
the conda `documentation <https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html>`_.

Package dependencies
--------------------

The Sockeye pipeline makes use of the following dependencies. No manual
installation is required, as these are all installed automatically into a series
of ``conda`` environments that are created throughout the course of a pipeline
run:

- bedtools [1_]
- bioframe [2_]
- biopython [3_]
- editdistance [4_]
- matplotlib [5_]
- minimap2 [6_]
- numpy [7_]
- pandas [8_]
- parasail-python [9_]
- pysam [10_]
- samtools [11_]
- scikit-learn [12_]
- seqkit [13_]
- tqdm [14_]
- umap-learn [15_]
- vsearch [16_]

Additionally, while no explicit dependency exists for the
`UMI-tools <https://github.com/CGATOxford/UMI-tools>`_ package  [17_], the script
``bin/cluster_umis.py`` makes significant use of several functions from
the package. More detailed acknowledgements can be found in the source code.
