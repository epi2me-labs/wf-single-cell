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

process call_adapter_scan {
    // Stranding of reads and ?
    // Neil: Only one thread for this. Seems low
    // Do we need a batch number to add to output filenames?
    input:
        tuple val(sample_id), path(fastq_chunk)
        tuple val(sample_id), val(kit)
    output:
        tuple val(sample_id), path("*.fastq"), emit: STRANDED_FQ_CHUNKED
        tuple val(sample_id), path("*.tsv"), emit: READ_CONFIG_CHUNKED
    """
    python adapter_scan_vsearch.py \
    -t 1 \
    --kit $kit \
    --output_fastq ${sample_id}_ad_scan.fastq \
    --output_tsv  ${sample_id}_ad_scan.tsv \
    --batch_size $params.READ_STRUCTURE_BATCH_SIZE \
    $fastq_chunk
    """
}



// process combine_adapter_tables: Do this with workflow operators
//combine_stranded_fastqs: Do this woth a workflow operator

process summarize_adapter_table {
    input:
        tupe val(sample_id), path(READ_CONFIG)
    output:
        tuple val(sample_id), path('*config_stats.json'), emit: CONFIG_STATS
    """
    #!/usr/bin/env python
    import pandas as pd
    import json

    df = pd.read_csv(READ_CONFIG, sep="\t")
    stats = {}
    stats[params.run_id] = {}
    stats[params.run_id]["general"] = {}
    stats[params.run_id]["general"]["n_reads"] = df.shape[0]
    stats[params.run_id]["general"]["rl_mean"] = df["readlen"].mean()
    stats[params.run_id]["general"]["rl_std_dev"] = df["readlen"].std()
    stats[params.run_id]["general"]["n_fl"] = df[df["fl"] == True].shape[0]
    stats[params.run_id]["general"]["n_stranded"] = df[
        df["stranded"] == True
    ].shape[0]

    stats[params.run_id]["strand_counts"] = {}
    stats[params.run_id]["strand_counts"]["n_plus"] = df[
        df["orig_strand"] == "+"
    ].shape[0]
    stats[params.run_id]["strand_counts"]["n_minus"] = df[
        df["orig_strand"] == "-"
    ].shape[0]

    stats[params.run_id]["detailed_config"] = {}
    for category, n in df["orig_adapter_config"].value_counts().items():
        stats[params.run_id]["detailed_config"][category] = n

    stats[params.run_id]["summary_config"] = {}
    for label, n in df["lab"].value_counts().items():
        stats[params.run_id]["summary_config"][label] = n

    with open(${sample_id}_config_stats.json, "w") as f:
        json.dump(stats, f, indent=4)
    """

}


// workflow module
workflow stranding {
    take:
        fastq
        refBase
        genome
        annotation
    main:
        chunk_files(fastq)
        call_adapter_scan(chunk_files.out)
        STRANDED_FQ = call_adapter_scan.out.STRANDED_FQ_CHUNKED.collectFile(keepHeader:false)
        summarize_adapter_table(
            call_adapter_scan.out.READ_CONFIG_CHUNKED.collectFile(keepHeader:true))
    emit:
        STRANDED_FQ = STRANDED_FQ
        CONFIG_STATS = summarize_adapter_table.out.CONFIG_STATS
}
