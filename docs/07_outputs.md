Outputs files may be aggregated including information for all             samples or provided per sample. Per sample files             will be prefixed with respective aliases and represented             below as {{ alias }}.

| Title | File path | Description | Per sample or aggregated |
|-------|-----------|-------------|--------------------------|
| workflow report | ./wf-single-cell-report.html | Report for all samples | aggregated |
| Concatenated sequence data | ./fastq_ingress_results/reads/{{ alias }}.fastq.gz | Per sample reads concatenated in a single FASTQ file. | per-sample |
| Results summaries | ./{{ alias }}/{{ alias }}.config_stats.json | Results summaries including adapter configuration numbers. | per-sample |
| Gene expression counts | ./{{ alias }}/{{ alias }}.gene_expression.counts.tsv | Gene x cell expression matrix. | per-sample |
| Processed gene expression counts | ./{{ alias }}/{{ alias }}.gene_expression.processed.tsv | Filtered and normalized gene x cell expression matrix. | per-sample |
| Transcript expression counts | ./{{ alias }}/{{ alias }}.transcript_expression.counts.tsv | Transcript x cell expression matrix. | per-sample |
| Processed transcript expression counts | ./{{ alias }}/{{ alias }}.transcript_expression.processed.tsv | Filtered and normalized transcript x cell expression matrix. | per-sample |
| Mitochondrial expression levels | ./{{ alias }}/{{ alias }}.gene_expression.mito.tsv | Per cell mitochondrial gene expression as percentage total of total gene expression. | per-sample |
| Knee plot | ./{{ alias }}/{{ alias }}.kneeplot.png | Knee plot illustrating the filtering of cells by read count. | per-sample |
| Saturation curves | ./{{ alias }}/{{ alias }}.saturation_curves.png | Saturation plots that indicate sampling of library complexity. | per-sample |
| Read tags | ./{{ alias }}/{{ alias }}.read_tags.tsv | Per read assigned barcodes UMIs genes and transcripts. | per-sample |
| Barcode counts | ./{{ alias }}/{{ alias }}.uncorrected_bc_counts.tsv | The counts of each barcode present in the sequenced library (only barcodes that have a 100% match in the 10x whitelist are included). | per-sample |
| Whitelist | ./{{ alias }}/{{ alias }}.whitelist.tsv | The barcodes found in the library that remain after filtering. | per-sample |
| Alignment output per chromosome | ./{{ alias }}/bams/{{ alias }}.{{ chromosome }}.tagged.bam | Genomic alignment output file per chromosome. | per-sample |
| Alignment index per chromosome | ./{{ alias }}/bams/{{ alias }}.{{ chromosome }}.tagged.bam.bai | Genomic alignment index file per chromosome. | per-sample |
| Alignment output per sample | ./{{ alias }}/bams/{{ alias }}.tagged.sorted.bam | Genomic alignment output file with aggregated chromosomes (when using --merge_bam). | per-sample |
| Alignment index per sample | ./{{ alias }}/bams/{{ alias }}.tagged.sorted.bam.bai | Genomic alignment index file with aggregated chromosomes (when using --merge_bam). | per-sample |