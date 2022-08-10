


process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "wftemplate"
    cpus params.threads
    input:
        tuple val(sample_id), 
              path(bam_sort), 
        val kit
        val barcode_length 
        val umi_length 
        path barcode_superlist
    output:
        tuple val(sample_id), path("*output.bam"), emit: bam
        tuple val(sample_id), path("*output.tsv"), emit: counts
    """
    samtools index $bam_sort
        extract_barcode.py \
        -t $task.cpus \
        --kit $kit \
        --adapter1_suff_length $params.BARCODE_ADAPTER1_SUFF_LENGTH \
        --barcode_length $barcode_length \
        --umi_length $umi_length \
        --output_bam "${sample_id}.output.bam" \
        --output_barcodes "${sample_id}.output.tsv" \
        $bam_sort $barcode_superlist
    """
}

process cleanup_headers_1{
    label "wftemplate"
    cpus params.threads
    input:
        tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bai"), emit: rehead
    """
        samtools reheader --no-PG -c 'grep -v ^@PG' $bam > "$sample_id".bam
        samtools index "$sample_id".bam
    
    """
}

process generate_whitelist{
    label "wftemplate"
    cpus params.threads
    input:
        tuple val(sample_id), path(counts)
    output:
        tuple val(sample_id), path("*whitelist.tsv"), emit: whitelist
        tuple val(sample_id), path("*kneeplot.png"), emit: kneeplot
    """
    knee_plot.py $params.barcode_kneeplot_flags --output_whitelist "$sample_id"_whitelist.tsv \
    --output_plot "$sample_id"_kneeplot.png $counts
    """
}

process split_bam_by_chroms{
    label "wftemplate"
    cpus params.threads
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        path "$sample_id/chr*.bam", emit: bams
        path "$sample_id/chr*.bai", emit: bai
    """
        split_bam_by_chroms.py -t ${task.cpus} --output_dir $sample_id $bam
    """
}



process assign_barcodes{
    label "wftemplate"
    cpus 10
    input:
        tuple val(sample_id), val(chrom), path(bam), path(bai)
        tuple val(sample_id), path(whitelist)
    output:
    """
    assign_barcodes.py -t ${task.cpus} --output_bam output --output_counts counts \
    --max_ed $params.max_ed --min_ed_diff $params.min_ed_diff --kit "$params.kit" \
    --adapter1_suff_length $params.adapter1_suff_length --barcode_length $params.barcode_length \
    --umi_length $params.umi_length $bam $whitelist
"""
}


process get_kit_info {
    input:
        path kit_config
        path sc_sample_sheet
              
    output:
        path 'sample_kit_info.csv'
    
    script:
    def bc_longlist_dir = "${projectDir}/data"
    """
    #!/usr/bin/env python
    import pandas as pd

    sample_df = pd.read_csv("${sc_sample_sheet}", sep=",", comment="#").set_index("run_id", drop=True)

    kit_df = pd.read_csv("$kit_config", sep=",", comment="#")

    records = []
    for run_id, row in sample_df.iterrows():
        kit_name = sample_df.loc[run_id, "kit_name"]
        kit_version = sample_df.loc[run_id, "kit_version"]

        # Get barcode length based on the kit_name and kit_version specified for this run_id.
        rows = (kit_df["kit_name"] == kit_name) & (kit_df["kit_version"] == kit_version)
        barcode_length = kit_df.loc[rows, "barcode_length"].values[0]

        # Get the appropriate cell barcode longlist based on the kit_name specified for this run_id.
        if kit_name == "3prime":
            long_list = "${bc_longlist_dir}/3M-february-2018.txt.gz"
        elif kit_name == "5prime":
            long_list = "${bc_longlist_dir}/737K-august-2016.txt.gz"
        elif kit_name == "multiome":
            long_list = "${bc_longlist_dir}/737K-arc-v1.txt.gz"
        else:
            raise Exception("Encountered an unexpected kit_name in samples.csv")

        # Get UMI length based on the kit_name and kit_version specified for this run_id.
        umi_length = kit_df.loc[rows, "umi_length"].values[0]
        records.append([run_id, kit_name, kit_version, barcode_length, umi_length, long_list])
    df_out = pd.DataFrame.from_records(records)
    df_out.columns = ["sample_id", "kit_name", "kit_version", "barcode_length", "umi_length", "long_list"]
    df_out.to_csv('sample_kit_info.csv', index=False)
    """
}

workflow process_bams {
    take:
        bam
        bam_idx
        sc_sample_sheet
        kit_config
    main:
        // bam_list = file("$projectDir/data/data/3M-february-2018.txt.gz")
        // bam_sort = file(params.bam_sort)
        // kit = params.kit
        // barcode_length = params.barcode_length
        // umi_length = params.umi_length
        // sample_id = "test"

        kit_details = get_kit_info(
            kit_config,
            sc_sample_sheet)
        
        // bam.join(kit_details).view()

        // extract_barcodes(
        //     bam.join(kit_details))
        
        // headers = cleanup_headers_1(newl.bam)
        // white_list = generate_whitelist(newl.counts)
        // by_chrom = split_bam_by_chroms(headers.rehead)
        // bam = by_chrom.bams.flatten().map{ 
        //     it -> tuple(it.toString().split('/')[-2].toString(), it.toString().split("/")[-1].split("\\.sorted")[0], it)}
        // bai = by_chrom.bai.flatten().map{ 
        //     it -> tuple(it.toString().split('/')[-2].toString(), it.toString().split("/")[-1].split("\\.sorted")[0], it)}
        // bam_bai = bam.join(bai, by:[0,1])
        // assign_barcodes(bam_bai,white_list.whitelist)
}