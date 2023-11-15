# wf-single-cell: Limitations and known issues  


## No trimmed FASTQ output 
Users may want to obtained FASTQ files trimmed of adapter sequences, barcodes and UMIs for their
downstream analysis, but there is not currently an option to output such files.

## No ability to filter non-full-length reads
Subreads are classified as being full length if flanked by two compatible adapters. 
However, at the moment this classification has no effect 
on whether these are further processed by the workflow. A user option to control this behaviour may be desirable. 

## 10x gene expression and feature barcodes discrimination
A 10x barcode whitelist containing all possible barcodes is used to cross-reference the discovered baroccodes for
barcode error correction. For the 3prime and multiome kits, the whitelist contains ~3M gene expression barcodes that this workflow is interested in.
However, it also contains a similar number of feature barcodes, see this [10x article](https://kb.10xgenomics.com/hc/en-us/articles/360031133451-Why-is-there-a-discrepancy-in-the-3M-february-2018-txt-barcode-whitelist-).
It's not currently possible to differentiate between the two types of barcode in this whitelist. 
Therefore, it is possible that somen error-containing gene expression barcodes may be being incorrectly assigned
to feature barcodes. To what extent this is happening is currently unknown. 

## Gene and feature assignment discrepancy 
In the `assign_features` process, genes are only assigned if they have a MAPQ score greater than a user-defined MAPQ score (default 30)
However,  transcripts are assigned based on alignment to a transcriptome that is built during the workflow. 
Transcripts are not filtered by MAPQ, but by applying some alternative heuristics based on alignment scores as well as transcript and query 
coverages. This can lead to cases where transcripts are called, but not genes. This will be fixed in a future version.
  