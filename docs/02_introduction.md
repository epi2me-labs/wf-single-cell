This workflow extracts cell barcodes and UMIs from [10x](https://www.10xgenomics.com/)-generated single cell libraries.
It was initially created as a Nextflow port of [Sockeye](https://github.com/nanoporetech/sockeye).

In brief, the workflow does the following:

+ Adapter identification, fused read splitting and stranding.
+ Mapping of reads to genomic reference.
+ Gene and transcript read assignment.
+ Cell barcode and UMI extraction and correction.
+ Generation of gene and transcript count matrices for unique UMIs.
+ Tagging BAM files with cell barcodes and UMIs.
+ Calculation of library saturation.

This workflow supports the following 10x kits:
+ 3': v2/v3 and v4 (GEM-X)
+ 5': v1/v2
+ multiome (gene expression only): v1 
+ visium spatial transcriptomics: v1


The [BLAZE](https://github.com/shimlab/BLAZE) preprint provided useful benchmarking of the original sockeye implementation. 
This assisted in the selection of appropriate thresholds for cell cut-off and for defining the limits of the gene x cell matrix.

The isoform selection procedure used in this workflow was adapted from that found in the [FLAMES](https://github.com/LuyiTian/FLAMES) 
package.
