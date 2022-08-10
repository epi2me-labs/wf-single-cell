process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "wftemplate"
    cpus params.threads
    input:
        tuple val(sample_id), path(bam_sort), val(kit), val(barcode_length), val(umi_length), path(list)
    output:
        tuple val(sample_id), path("*output.bam"), emit: bam
        tuple val(sample_id), path("*output.tsv"), emit: counts
    """
    samtools index $bam_sort
    extract_barcode.py -t $task.cpus --kit $kit --adapter1_suff_length $params.adapter1_suff_length \
    --barcode_length $barcode_length --umi_length $umi_length --output_bam "$sample_id".output.bam \
    --output_barcodes "$sample_id".output.tsv $bam_sort $list
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

workflow process_bams {
    take:
        bam
    main:
        bam_list = file("$projectDir/data/data/3M-february-2018.txt.gz")
        bam_sort = file(params.bam_sort)
        kit = params.kit
        barcode_length = params.barcode_length
        umi_length = params.umi_length
        sample_id = "test"
        newl = extract_barcodes([sample_id, bam_sort, kit, barcode_length, umi_length,bam_list])
        headers = cleanup_headers_1(newl.bam)
        white_list = generate_whitelist(newl.counts)
        by_chrom = split_bam_by_chroms(headers.rehead)
        bam = by_chrom.bams.flatten().map{ 
            it -> tuple(it.toString().split('/')[-2].toString(), it.toString().split("/")[-1].split("\\.sorted")[0], it)}
        //bams.view()
        bai = by_chrom.bai.flatten().map{ 
            it -> tuple(it.toString().split('/')[-2].toString(), it.toString().split("/")[-1].split("\\.sorted")[0], it)}
        //b//am.view()
        //bai.view()
        lol = bam.join(bai, by:[0,1])
        assign_barcodes(lol,white_list.whitelist)

    emit:
      kit
}