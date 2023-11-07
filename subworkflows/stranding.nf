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
        tuple val(meta),
              path(chunk, stageAs: 'chunk.fq.gz')
    output:
        tuple val(meta), path("${meta.alias}_adapt_scan.fastq.gz"), emit: stranded_fq_chunked
        tuple val(meta), path("${meta.alias}_adapt_scan.tsv"), emit: read_config_chunked
    
    """
    export POLARS_MAX_THREADS=$task.cpus

    workflow-glue adapter_scan_vsearch \
    chunk.fq.gz \
    --kit ${meta['kit_name']} \
    --output_fastq "${meta.alias}_adapt_scan.fastq.gz" \
    --output_tsv  "${meta.alias}_adapt_scan.tsv"
    """
}


process combine_adapter_tables {
    label "singlecell"
    cpus 1
    input:
        tuple val(meta), path("adapters.tsv")
    output:
        tuple val(meta), path("${meta.alias}_read_config.tsv"), emit: read_config
    """
    # Concatenate tsv file keeping header from first file.
    awk 'FNR==1 && NR!=1{next;}{print}' adapters.tsv* > "${meta.alias}_read_config.tsv"
    """
}


process summarize_adapter_table {
    label "singlecell"
    cpus 1
    memory 1.5.GB
    input:
        tuple val(meta), path(read_config)
    output:
        tuple val(meta), path("${meta.alias}.config_stats.json"), emit: config_stats
    """
    workflow-glue summarise_adapters \
        --read_config_tsv "${read_config}" \
        --sample_id "${meta.alias}" \
        --out "${meta.alias}.config_stats.json" \
        --threads $task.cpus
    """
}


// workflow module
workflow stranding {
    take:
        read_chunks
    main:
         // Rejig checks to so each is  [meta, fastq_chunk]
        meta_chunks = read_chunks.flatMap({it ->
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        })
        call_adapter_scan(meta_chunks)
        combine_adapter_tables(call_adapter_scan.out.read_config_chunked.groupTuple())
        summarize_adapter_table(combine_adapter_tables.out.read_config)

    emit:
        stranded_fq = call_adapter_scan.out.stranded_fq_chunked
        config_stats = summarize_adapter_table.out.config_stats
}
