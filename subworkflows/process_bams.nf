


process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "wftemplate"
    cpus params.max_threads
    input:
        tuple val(sample_id), 
              val(kit_name),
              val(kit_version),
              val(barcode_length),
              val(umi_length),
              val(bc_long_list),
              path(bam),
              path(bam_idx)

    output:
        tuple val(sample_id), path("*.bam"), path("*.bai"), emit: bam
        tuple val(sample_id), path("*.tsv"), emit: counts
    """
    extract_barcode.py \
    $bam $bc_long_list\
    -t $task.cpus \
    --kit $kit_name \
    --adapter1_suff_length $params.BARCODE_ADAPTER1_SUFF_LENGTH \
    --barcode_length $barcode_length \
    --umi_length $umi_length \
    --output_bam "tmp.bam" \
    --output_barcodes "${sample_id}.barcodes.tsv";

    samtools reheader --no-PG -c 'grep -v ^@PG' tmp.bam > "${sample_id}.bam"
    samtools index "${sample_id}.bam"
    rm tmp.bam
    """
}

process generate_whitelist{
    label "wftemplate"
    cpus params.max_threads
    input:
        tuple val(sample_id), path(counts)
    output:
        tuple val(sample_id), path("*whitelist.tsv"), emit: whitelist
        tuple val(sample_id), path("*kneeplot.png"), emit: kneeplot
    """
    knee_plot.py ${params.BARCODE_KNEEPLOT_FLAGS} \
        --output_whitelist "${sample_id}_whitelist.tsv" \
        --output_plot "${sample_id}_kneeplot.png" $counts
    """
}

process split_bam_by_chroms{
    label "wftemplate"
    cpus params.max_threads
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        tuple val(sample_id), 
              path("splits/chr*.bam"),
              path("splits/chr*.bai"), emit: bam
    """
    split_bam_by_chroms.py -t ${task.cpus} --output_dir splits $bam
    """
}


process assign_barcodes{
    label "wftemplate"
    cpus params.max_threads
    input:
         tuple val(sample_id), 
               val(kit_name),
               val(kit_version),
               val(barcode_length),
               val(umi_length),
               path(bc_long_list),
               path(whitelist),
               val(chr),
               val(sample_id),
               path(bam),
               path(bai)
    output:
        tuple val(sample_id), 
            val(chr),
            path("${sample_id}.bc_assign.bam"),
            path("${sample_id}.bc_assign.bam.bai"), 
            emit: CHROM_BAM_BC_BAI
        tuple val(sample_id), 
              path("${sample_id}.bc_assign_counts.tsv"), 
              emit: CHROM_ASSIGNED_BARCODE_COUNTS
    """
    assign_barcodes.py -t ${task.cpus} \
        --output_bam tmp.bam \
        --output_counts ${sample_id}.bc_assign_counts.tsv \
        --max_ed $params.BARCODE_MAX_ED \
        --min_ed_diff $params.BARCODE_MIN_ED_DIFF \
        --kit $kit_name \
        --adapter1_suff_length $params.BARCODE_ADAPTER1_SUFF_LENGTH \
        --barcode_length $barcode_length \
        --umi_length \
        $umi_length $bam $whitelist

    # Cleanup header #2
    samtools reheader --no-PG -c 'grep -v ^@PG' \
        tmp.bam > "${sample_id}.bc_assign.bam";
    samtools index "${sample_id}.bc_assign.bam"
"""
}

process bam_to_bed {
    input:
        tuple val(sample_id), 
              val(chr),
              path(bam),
              path(bai) 
    output:
        tuple val(sample_id), 
              val(chr),
              path('*bc_assign.bed'), 
              emit: CHROM_BED_BC
    """
    bedtools bamtobed -i $bam > "${sample_id}_bc_assign.bed"
    """
}

process split_gtf_by_chroms {
    input:
        path(gtf)
    output:
        path("*"), emit: CHROM_GTF
    """
    awk '/^[^#]/ {print>\$1}' $gtf 
    """
}   

process assign_genes {
    input:
        tuple val(sample_id),
              val(chr),
              path(CHROM_BED_BC),
              path(CHROM_GTF)
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_${chr}.read.gene_assigns.tsv"),
              emit: CHROM_TSV_GENE_ASSIGNS
    """
    python assign_genes.py \
    --output ${sample_id}_${chr}.read.gene_assigns.tsv \
    CHROM_BED_BC CHROM_GTF
    """
}

process add_gene_tags_to_bam {
    input:
        tuple val(sample_id),
              cal(chr),
              path(CHROM_BAM_BC),
              path(CHROM_BAM_BC_BAI),
              path(CHROM_TSV_GENE_ASSIGNS)
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_bc_assign.gene.bam"), 
              path("${sample_id}_bc_assign.gene.bam"),
              emit: CHROM_BAM_BC_GENE_BAI
    """
    touch {input.bai};
    add_gene_tags.py \
        --output tmp.bam \
    CHROM_BAM_BC CHROM_TSV_GENE_ASSIGNS"

    # Cleanup headers #3
    samtools reheader --no-PG -c 'grep -v ^@PG' \
    tmp.bam > ; ${sample_id}_bc_assign.gene.bam
    samtools index ${sample_id}_bc_assign.gene.bam
    """
}


process cluster_umis {
    cpus params.UMI_CLUSTER_MAX_THREADS
    input:
        tuple val(sample_id),
              val(chr),
              path(bam),
              path(bai),
    output:
         tuple val(sample_id),
              val(chr),
              path("*.bam"),
              path("*.bam.bai")
    // params:
    //     interval=config["UMI_GENOMIC_INTERVAL"],
    //     cell_gene_max_reads=config["UMI_CELL_GENE_MAX_READS"],
    // threads: config["UMI_CLUSTER_MAX_THREADS"]
    """
    cluster_umis.py 
    --threads task.cpus \
    --ref_interval $params.UMI_GENOMIC_INTERVAL \ 
    --cell_gene_max_reads $params.UMI_CELL_GENE_MAX_READS \
    --output tmp.bam $bam

    //Cleanup headers #4
    samtools reheader --no-PG -c 'grep -v ^@PG' \
    tmp.bam > {output.bam};
    samtools index ${sample_id}_${chr}_tagged.bam
    """
}


process get_kit_info {
    input:
        path kit_config
        path sc_sample_sheet
              
    output:
        path 'sample_kit_info.csv', emit: kit_info
    
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
        gtf
    main:
        get_kit_info(
            kit_config,
            sc_sample_sheet)


        extract_barcodes(
            get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(bam).join(bam_idx))
        
        split_bam_by_chroms(
            extract_barcodes.out.bam
        )

        generate_whitelist(
            extract_barcodes.out.counts
        )

        bam_bai_chromes = split_bam_by_chroms.out.bam.map({it ->
        // Merge bams and indxes along and add chromosome into tuple
        // [sample_id, chr, bam, bai]
            pairs = []
            for (i=0; i<it[1].size(); i++) {
                chr = it[1][i].toString().tokenize('/')[-1].tokenize('.')[-3]
                pairs.add(tuple(it[0], chr, it[1][i], it[2][i]))
            }
            return pairs
        }).flatMap(it-> it)

        // // Todo remove redundant sample_ID
        chro_bam_kit = get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1).join(generate_whitelist.out.whitelist)
            .cross(bam_bai_chromes).map({it -> it.flatten()})
            
        assign_barcodes(chro_bam_kit)

        bam_to_bed(assign_barcodes.out.CHROM_BAM_BC_BAI)

        chr_gtf = split_gtf_by_chroms(gtf)
        .flatten()
        .map {file -> 
            // create [chr, gtf]
             tuple(file.toString().tokenize('/')[-1], file)}.view() // Add chromosome to tuple
           
 
        // assign_genes(
        //     bam_to_bed.out.CHROM_BED_BC
        //     .cross(chr_gtf{it -> it[1]})
        //     )


        
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