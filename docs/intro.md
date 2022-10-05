## Introduction

The following single-cell kits from 10x Genomics are currently supported:
- Chromium Single Cell [3ʹ gene expression](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium3.html), versions 2 and 3
- Chromium Single Cell [5ʹ gene expression](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium5.html), version 1
- Chromium Single Cell [Multiome (ATAC + GEX)](https://teichlab.github.io/scg_lib_structs/methods_html/10xChromium_multiome.html), version 1

Oxford Nanopore has developed a protocol for sequencing single-cell libraries from 10x, which can be found on the Nanopore Community [website](https://community.nanoporetech.com/docs/prepare/library_prep_protocols/single-cell-transcriptomics-10x/v/sst_v9148_v111_revb_12jan2022).

The inputs to Sockeye are raw nanopore reads (FASTQ) generated from the sequencing
instrument and reference files that can be downloaded from [10x](https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/latest).
The pipeline output a gene x cell, and transcript x cell expression matrices, as well as a BAM file of
aligned reads tagged with cell barcode and UMI information.

The BLAZE preprint provided useful benchmarking of the original sockeye implementation. This assisted in the selection of appropriate parameters for cell cut-off thresholds and for defining the limits of the cell x gene matrix.