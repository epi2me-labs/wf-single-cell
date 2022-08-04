// This rule looks for all files in a input directory. No need to do that as fastcat
// has already concatenated all reads
//rule cp_batch_fastqs:


process chunk_files {
    // The orginal SM rule (call_cat_fastq) called chunk_fastqs.py, which did file concatenation and chuck creation, we need only the latter
    input:
        tuple val(sample_id), path(fastq)
    output:
        tuple val(sample_id), path("*.fastq")
    """
    seqkit split $fastq -p --two-pass $params.MAX_THREADS
    """
}

process adapter_scan {
    // Neil: Only one thread for this. Seems low
    // Do we need a batch number to add to output filenames?
    input:
        tuple val(sample_id), path(fastq_chunk)
        tuple val(sample_id), val(kit)
    output:
        tuple val(sample_id), path("*.fastq"), emit: fastq
        tuple val(sample_id), path("*.tsv"), emit: tsv
    """
    python adapter_scan_vsearch.py \
    -t 1 \
    --kit $kit \
    --output_fastq ${sample_id}_ad_scan.fastq \
    --output_tsv  ${sample_id}_add_scan.tsv \
    --batch_size $params.READ_STRUCTURE_BATCH_SIZE \
    $fastq_chunk
    """
}