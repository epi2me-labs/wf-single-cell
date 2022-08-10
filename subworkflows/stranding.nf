// This rule looks for all files in a input directory. No need to do that as fastcat
// has already concatenated all reads
//rule cp_batch_fastqs:


process chunk_files {
    // The orginal SM rule (call_cat_fastq) called chunk_fastqs.py, which did file concatenation and chuck creation, we need only the latter
    input:
        tuple val(sample_id),
              val(kit_name),
              val(kit_version), 
              path(fastq)

    output:
        tuple val(sample_id),
              val(kit_name),
              val(kit_version),
              path("chunks/*")
    """
    seqkit split $fastq -p $params.MAX_THREADS -O chunks
    """
}

process call_adapter_scan {
    // Stranding of reads and ?
    // Neil: Only one thread for this. Seems low
    // Do we need a batch number to add to output filenames?
    input:
        tuple val(sample_id),
            val(kit_name),
            val(kit_version), 
            path(fastq_chunk)
    output:
        tuple val(sample_id), path("*.fastq"), emit: STRANDED_FQ_CHUNKED
        tuple val(sample_id), path("*.tsv"), emit: READ_CONFIG_CHUNKED
    """
    # Get batch name from file to prevent name collisions
    batch=\$(echo ${fastq_chunk}|awk -F'.' '{print \$(NF-2)}')
    echo \$batch
    
    adapter_scan_vsearch.py \
    $fastq_chunk \
    -t 1 \
    --kit $kit_name \
    --output_fastq ${sample_id}_\${batch}_adapt_scan.fastq \
    --output_tsv  ${sample_id}_\${batch}_adapt_scan.tsv \
    --batch_size $params.READ_STRUCTURE_BATCH_SIZE \
    """
}

// process combine_adapter_tables: Do this with workflow operators
//combine_stranded_fastqs: Do this woth a workflow operator

process combine_adapter_tables {
    input:
        tuple val(sample_id), path(tsv_files)
    output:
        tuple val(sample_id), path("*read_config.tsv"), emit: READ_CONFIG
    """
    # Concatenate tsv file keeping header from first file.
    awk 'FNR==1 && NR!=1{next;}{print}' *.tsv > ${sample_id}_read_config.tsv
    """
}

process gather_fastq{
    input:
        tuple val(sample_id), path(fastq_files)
    output:
        tuple val(sample_id), path('*merged.fastq'), emit: merged_fastq
    """
    cat *.fastq >> ${sample_id}_merged.fastq
    """
}

process summarize_adapter_table {
    input:
        tuple val(sample_id), path(READ_CONFIG)
    output:
        tuple val(sample_id), path('*config_stats.json'), emit: CONFIG_STATS
    """
    #!/usr/bin/env python
    import pandas as pd
    import json

    df = pd.read_csv("${READ_CONFIG}", sep="\t")
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

    with open("${sample_id}_config_stats.json", "w") as f:
        json.dump(stats, f, indent=4)
    """

}


// workflow module
workflow stranding {
    take:
        inputs
        sc_sample_sheet
    main:
        d = {it ->
        /* Harmonize tuples
        output:
            tuple val(sample_id), path('*.gff')
        When there are multiple paths, will emit:
            [sample_id, [path, path ..]]
        when there's a single path, this:
            [sample_id, path]
        This closure makes both cases:
            [[sample_id, path][sample_id, path]].
        */
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }

        chunk_files(inputs)
        
        chunk_files.out.flatMap(d)
        
        call_adapter_scan(
            chunk_files.out.flatMap({
                // TODO: look for a better way to emit files in tuples using wildcards
                if (it[3].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
                }
                l = [];
                for (x in it[3]){
                    l.add(tuple(it[0], it[1], it[2], x))
                }
                return l
            }))
        
        gather_fastq(call_adapter_scan.out.STRANDED_FQ_CHUNKED.groupTuple())

        combine_adapter_tables(call_adapter_scan.out.READ_CONFIG_CHUNKED.groupTuple())
        summarize_adapter_table(combine_adapter_tables.out.READ_CONFIG)
            
    emit:
        STRANDED_FQ = gather_fastq.out.merged_fastq
        CONFIG_STATS = summarize_adapter_table.out.CONFIG_STATS
}
