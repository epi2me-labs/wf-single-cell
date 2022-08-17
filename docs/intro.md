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
