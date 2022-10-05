process call_adapter_scan {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              val(meta),
              path(fastq_chunk)
    output:
        tuple val(sample_id), path("*.fastq"), emit: stranded_fq_chunked
        tuple val(sample_id), path("*.tsv"), emit: read_config_chunked
    
    """
    # Get batch name from file to prevent name collisions
    batch=\$(echo ${fastq_chunk}|awk -F'.' '{print \$(NF-2)}')
    echo \$batch
    
    adapter_scan_vsearch.py \
    $fastq_chunk \
    -t 1 \
    --kit ${meta['kit_name']} \
    --output_fastq ${sample_id}_\${batch}_adapt_scan.fastq \
    --output_tsv  ${sample_id}_\${batch}_adapt_scan.tsv \
    --batch_size $params.read_structure_batch_size \
    """
}

// process combine_adapter_tables: Do this with workflow operators
//combine_stranded_fastqs: Do this woth a workflow operator

process combine_adapter_tables {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), path(tsv_files)
    output:
        tuple val(sample_id), path("*read_config.tsv"), emit: read_config
    """
    # Concatenate tsv file keeping header from first file.
    awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > ${sample_id}_read_config.tsv
    """
}

process gather_fastq{
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), path(fastq_files)
    output:
        tuple val(sample_id), path('*merged.fastq'), emit: merged_fastq
    """
    cat *.fastq >> ${sample_id}_merged.fastq
    """
}

process summarize_adapter_table {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), path(read_config)
    output:
        tuple val(sample_id), path('*config_stats.json'), emit: config_stats
    """
    #!/usr/bin/env python
    import pandas as pd
    import json

    df = pd.read_csv("${read_config}", sep="\t")
    stats = {}
    stats["{$sample_id}"] = {}
    stats["{$sample_id}"]["general"] = {}
    stats["{$sample_id}"]["general"]["n_reads"] = df.shape[0]
    stats["{$sample_id}"]["general"]["rl_mean"] = df["readlen"].mean()
    stats["{$sample_id}"]["general"]["rl_std_dev"] = df["readlen"].std()
    stats["{$sample_id}"]["general"]["n_fl"] = df[df["fl"] == True].shape[0]
    stats["{$sample_id}"]["general"]["n_stranded"] = df[
        df["stranded"] == True
    ].shape[0]

    stats["{$sample_id}"]["strand_counts"] = {}
    stats["{$sample_id}"]["strand_counts"]["n_plus"] = df[
        df["orig_strand"] == "+"
    ].shape[0]
    stats["{$sample_id}"]["strand_counts"]["n_minus"] = df[
        df["orig_strand"] == "-"
    ].shape[0]

    stats["{$sample_id}"]["detailed_config"] = {}
    for category, n in df["orig_adapter_config"].value_counts().items():
        stats["{$sample_id}"]["detailed_config"][category] = n

    stats["{$sample_id}"]["summary_config"] = {}
    for label, n in df["lab"].value_counts().items():
        stats["{$sample_id}"]["summary_config"][label] = n

    with open("${sample_id}.config_stats.json", "w") as f:
        json.dump(stats, f, indent=4)
    """
}


// workflow module
workflow stranding {
    take:
        read_chunks
        meta
    
    main:      
        chunks = meta
            .cross(read_chunks.flatMap({it ->
            // Rejig the outputs to be [sample_id, fatq_chunk]
            // Then merge in kit info
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }))
        // sample_id, meta, fastq_chunk
        .map {it -> tuple(it[0][0], it[0][1], it[1][1])}

        call_adapter_scan(chunks)

        gather_fastq(call_adapter_scan.out.stranded_fq_chunked.groupTuple())

        combine_adapter_tables(call_adapter_scan.out.read_config_chunked.groupTuple())
        summarize_adapter_table(combine_adapter_tables.out.read_config)
            
    emit:
        stranded_fq = gather_fastq.out.merged_fastq
        config_stats = summarize_adapter_table.out.config_stats
}
