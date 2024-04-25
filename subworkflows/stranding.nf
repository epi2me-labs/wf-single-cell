process call_adapter_scan {
    label "singlecell"
    cpus 2
    // Benchmarking has shown that memory usage is ~ 1.5x times fastq size.
    // Smaller chunks sizes have a larger ratios, so 1G is added to account for this.
    // Occasionally memory requirements are higher so attempt retries with increasing memory too.
    memory { 1.0.GB.toBytes() + (chunk.size() * 2) * task.attempt}
    maxRetries = 3
    errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        tuple val(group_index),
              val(meta),
              path(chunk, stageAs: 'chunk.fq.gz')
    output:
        tuple val(group_index),
              val(meta),
              path("adapt_scan.fastq.gz"),
              emit: stranded_fq_chunked
        tuple val(meta),
              path("adapt_scan.tsv"),
              emit: read_config_chunked
    script:
    def fl = params.full_length_only ? "--keep_fl_only": ""
    """
    export POLARS_MAX_THREADS=$task.cpus

    # writing fastq occurs in a single-loop - fine to give bgzip some more juice
    workflow-glue adapter_scan_vsearch \
        chunk.fq.gz \
        --kit ${meta['kit_name']} \
        --output_tsv  "adapt_scan.tsv" \
        ${fl} \
        | bgzip --threads ${task.cpus} -c > "adapt_scan.fastq.gz"
    """
}


process combine_adapter_tables {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta), path("adapters.tsv")
    output:
        tuple val(meta), path("read_config.tsv"), emit: read_config
    """
    # Concatenate tsv file keeping header from first file.
    awk 'FNR==1 && NR!=1{next;}{print}' adapters.tsv* > "read_config.tsv"
    """
}


process summarize_adapter_table {
    label "singlecell"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(read_config)
    output:
        tuple val(meta), path("config_stats.json"), emit: config_stats
    """
    workflow-glue summarise_adapters \
        --read_config_tsv "${read_config}" \
        --sample_id "${meta.alias}" \
        --out "config_stats.json" \
        --threads $task.cpus
    """
}


process extract_barcodes{
    label "singlecell"
    cpus 2
    memory "2.0 GB"
    input:
        tuple val(group_index),
              val(meta),
              path("fastq_chunk.fq")
        path "bc_longlist_dir"

    output:
        tuple val(group_index),
              val(meta),
              path("bc_extract.tsv"),
              emit: extracted_bc_umi
        tuple val(meta),
              path("high_quality_bc_counts.tsv"),
              emit: barcode_counts
        tuple val(group_index),
              val(meta),
              path("stranded_trimmed.fastq"),
              emit: trimmed_fq

    """
    workflow-glue extract_barcode \
        fastq_chunk.fq \
        bc_longlist_dir/${meta['bc_long_list']} \
        -t ${task.cpus} \
        --kit ${meta["kit_name"]} \
        --adapter1_suff_length $params.barcode_adapter1_suff_length \
        --min_barcode_qv $params.barcode_min_quality \
        --barcode_length ${meta['barcode_length']} \
        --umi_length ${meta['umi_length']} \
        --output_read_tags "bc_extract.tsv" \
        --output_barcode_counts "high_quality_bc_counts.tsv" \
        --output_trimmed_fastq "stranded_trimmed.fastq"
    """
}


// workflow module
workflow stranding {
    take:
        read_chunks
        bc_longlist_dir
    main:
        call_adapter_scan(read_chunks)
        combine_adapter_tables(
            call_adapter_scan.out.read_config_chunked.groupTuple())
        summarize_adapter_table(combine_adapter_tables.out.read_config)
        extract_barcodes(
            call_adapter_scan.out.stranded_fq_chunked,
            bc_longlist_dir)

    emit:
        stranded_trimmed_fq = extract_barcodes.out.trimmed_fq
        extracted_barcodes = extract_barcodes.out.extracted_bc_umi
        high_qual_bc_counts = extract_barcodes.out.barcode_counts
        config_stats = summarize_adapter_table.out.config_stats
}
