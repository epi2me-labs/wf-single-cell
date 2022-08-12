


process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "wfsockeye"
    conda "envs/barcodes.yml"
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
        tuple val(sample_id), path("*.bam"), emit: bam
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
    """
}

process cleanup_headers_1 {
    label "wfsockeye"
    cpus params.max_threads
    conda "envs/samtools.yml"
    input:
         tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bam.bai"), 
        emit: extract_barcodes_cleaned
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam > "${sample_id}.bam"
    samtools index "${sample_id}.bam"
    """
}

process generate_whitelist{
    label "wfsockeye"
    conda "envs/kneeplot.yml"
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
    label "wfsockeye"
    conda "envs/barcodes.yml"
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
    label "wfsockeye"
    conda "envs/barcodes.yml"
    cpus params.max_threads
    input:
         tuple val(sample_id), 
               val(kit_name),
               val(kit_version),
               val(barcode_length),
               val(umi_length),
               path(bc_long_list),
               path(whitelist),
               val(_), // Redundant sample_id: remove
               val(chr), 
               path(bam),
               path(bai)
    output:
        tuple val(sample_id), 
            val(chr),
            path("tmp.bam"),
            emit: CHROM_BAM_BC
        tuple val(sample_id),
              val(chr),
              path("*.bc_assign_counts.tsv"), 
              emit: CHROM_ASSIGNED_BARCODE_COUNTS
    """
    assign_barcodes.py -t ${task.cpus} \
        --output_bam tmp.bam \
        --output_counts ${sample_id}_${chr}.bc_assign_counts.tsv \
        --max_ed $params.BARCODE_MAX_ED \
        --min_ed_diff $params.BARCODE_MIN_ED_DIFF \
        --kit $kit_name \
        --adapter1_suff_length $params.BARCODE_ADAPTER1_SUFF_LENGTH \
        --barcode_length $barcode_length \
        --umi_length \
        $umi_length $bam $whitelist
"""
}

process cleanup_headers_2 {
    label "wfsockeye"
    cpus params.max_threads
    conda "envs/samtools.yml"
    input:
         tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), 
              val(chr), 
              path("*.bam"), 
              path("*.bam.bai"), 
        emit: CHROM_BAM_BC_BAI
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        > "${sample_id}.bc_assign.bam";
    samtools index "${sample_id}.bc_assign.bam"
    """
}

process bam_to_bed {
    label "wfsockeye"
    conda "envs/bedtools.yml"
    input:
        tuple val(chr), //emit chr first for doing cross on gtfs 
              val(sample_id),
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
    label "wfsockeye"
    conda "envs/assign_genes.yml"
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
    assign_genes.py \
    --output ${sample_id}_${chr}.read.gene_assigns.tsv \
    $CHROM_BED_BC $CHROM_GTF
    """
}

process add_gene_tags_to_bam {
    label "wfsockeye"
    conda "envs/umis.yml"
    input:
        tuple val(sample_id),
              cal(chr),
              path(CHROM_BAM_BC),
              path(CHROM_BAM_BC_BAI)
              path(CHROM_TSV_GENE_ASSIGNS)
    output:
        tuple val(sample_id),
              val(chr),
              path("tmp.bam"), 
              emit: CHROM_BAM_BC
    """
    add_gene_tags.py \
        --output tmp.bam \
        $CHROM_BAM_BC $CHROM_TSV_GENE_ASSIGNS
    """
}

process cleanup_headers_3 {
    label "wfsockeye"
    cpus params.max_threads
    conda "envs/samtools.yml"
    input:
         tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bam.bai"), 
        emit: CHROM_BAM_BC_BAI
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        > ${sample_id}_bc_assign.gene.bam;
    samtools index ${sample_id}_bc_assign.gene.bam
    """
}


process cluster_umis {
    label "wfsockeye"
    conda "envs/umis.yml"
    cpus params.UMI_CLUSTER_MAX_THREADS
    input:
        tuple val(sample_id),
              val(chr),
              path(bam),
              path(bai),
    output:
         tuple val(sample_id),
              val(chr),
              path("tmp.bam"),
              emit: bam
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
    """
}

process cleanup_headers_4{
    label "wfsockeye"
    conda "envs/samtools.yml"
    cpus params.max_threads
    input:
         tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bam.bai"), 
        emit: CHROM_BAM_BC_BAI
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        >  ${sample_id}_${chr}.tagged.bam;
    samtools index ${sample_id}_${chr}.tagged.bam
    """
}

process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "wfsockeye"
    conda "envs/samtools.yml"
    input:
        tuple val(sample_ids), path(bams)
    output:
        tuple val(sample_id), 
              path("${sample_id}.merged.bam"), 
              path("${sample_id}.merged.bam.bai"),
              emit: BAM_FULLY_TAGGED
    """
    samtools merge -o "${sample_id}.merged.bam" $bams; 
    samtools index "${sample_id}.merged.bam";
    """
}

process count_cell_gene_umi_reads {
    label "wfsockeye"
    conda "envs/barcodes.yml"
    input:
        tuple val(sample_id),
              path(bam)
    output:
        tuple val(sample_id)
              path("${sample_id}_cell_umi_gene.tsv"),
              emit: CELL_UMI_GENE_TSV
    """
    python cell_umi_gene_table.py \
        --output  ${sample_id}_cell_umi_gene.tsv $bam
    """
}

process umi_gene_saturation {
    label "wfsockeye"
    conda "envs/plotting.yml"
    input:
        tuple val(sample_id),
              path(cell_umi_gene_tsv)
    output:
        tuple val(sample_id),
              path("*saturation_curves.png")
    """
    python calc_saturation.py \
        --output ${sample_id}_saturation_curves.png" cell_umi_gene_tsv
    """
}

process construct_expression_matrix {
    label "wfsockeye"
    conda "envs/barcodes.yml"
    input:
        tuple val(sample_id),
              path(bam),
              path(bai)
    output:
        tuple val(sample_id), 
              path("*gene_expression.counts.tsv"), 
              emit: MATRIX_COUNTS_TSV
    """
    python gene_expression.py \
        --output ${sample_id}_gene_expression.counts.tsv $bam
    """
}

process process_expression_matrix {
    label "wfsockeye"
    conda "envs/barcodes.yml"
    input:
        tuple val(sample_id),
              path(MATRIX_COUNTS_TSV)
    output:
        tuple val(sample_id), 
        path("*gene_expression.processed.tsv"),
        emit: MATRIX_PROCESSED_TSV
    """
    process_matrix.py \
    --min_genes $params.MATRIX_MIN_GENES \
    --min_cells $params.MATRIX_MIN_CELLS \
    --max_mito $params.MATRIX_MAX_MITO \
    --norm_count $params.MATRIX_NORM_COUNT \
    --output ${sample_id}_gene_expression.processed.tsv $MATRIX_COUNTS_TSV
    """
}

process umap_reduce_expression_matrix {
    label "wfsockeye"
    conda "envs/umap.yml"
    input:
        tuple val(sample_id),
              path(MATRIX_PROCESSED_TSV)
    output:
         tuple val(sample_id),
              path("*gene_expression.umap.tsv"), 
              emit: MATRIX_UMAP_TSV
    """
    umap_reduce.py \
        --output ${sample_id}_gene_expression.umap.tsv $MATRIX_PROCESSED_TSV
    """
}


process umap_plot_total_umis {
    label "wfsockeye"
    conda "envs/plotting.yml"
    input:
        tuple val(sample_id),
              path(MATRIX_UMAP_TSV),
              path(MATRIX_PROCESSED_TSV)
    output:
          tuple val(sample_id),
              path("*umap.total.png"), 
              emit: MATRIX_UMAP_PLOT_TOTAL
    """
    python plot_umap.py \
        --output ${sample_id}_umap.total.png \
        $MATRIX_UMAP_TSV $MATRIX_PROCESSED_TSV
    """
}

process umap_plot_genes {
    // TODO: make a channle of input genes for thes process
    label "wfsockeye"
    conda "envs/plotting.yml"
    input:
        tuple val(sample_id),
              path(MATRIX_UMAP_TSV),
              path(MATRIX_PROCESSED_TSV),
              val(gene)
    output:
        tuple val(sample_id),
              path("*umap.gene.${gene}.png"), 
              emit: MATRIX_UMAP_PLOT_GENE
    """
    plot_umap.py \
        --gene $gene \
        --output {output.plot \
        $MATRIX_UMAP_TSV $MATRIX_PROCESSED_TSV
    """
}

process umap_plot_mito_genes {
    label "wfsockeye"
    conda "envs/plotting.yml"
    input:
        tuple val(sample_id),
              path(MATRIX_UMAP_TSV),
              path(MATRIX_PROCESSED_TSV)
    output:
        tuple val(sample_id),
              path("*umap.mitochondrial.png"), 
              emit: MATRIX_UMAP_PLOT_MITO
    conda "envs/plotting.yml"
    """
    plot_umap.py \
        --mito_genes \
        --output ${sample_id}_umap.mitochondrial.png \
        $MATRIX_UMAP_TSV $MATRIX_PROCESSED_TSV
    """
}
    

process get_kit_info {
    label "wfsockeye"
    conda "envs/assign_genes.yml"
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
        
        cleanup_headers_1(
            extract_barcodes.out.bam
        )
        
        split_bam_by_chroms(
            cleanup_headers_1.out.extract_barcodes_cleaned
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

        cleanup_headers_2(assign_barcodes.out.CHROM_BAM_BC)

        bam_to_bed(cleanup_headers_2.out.CHROM_BAM_BC_BAI)

        chr_gtf = split_gtf_by_chroms(gtf)
        .flatten()
        .map {file -> 
            // create [chr, gtf]
             tuple(file.toString().tokenize('/')[-1], file)} // Add chromosome to tuple

        assign_genes(
            bam_to_bed.out.CHROM_BED_BC
            .cross(chr_gtf).map({it ->
            // rejig the tuple to be [sample_id, chr, bed, gtf]
            [it[0][1], it[0][0], it[0][2], it[1][1]]})
        )

        add_gene_tags_to_bam(
            cleanup_headers_2.out.CHROM_BAM_BC_BAI
            // join on sample_id + chr
            .join(assign_genes.out.CHROM_TSV_GENE_ASSIGNS, by:[0, 1]))

        // cleanup_headers_3(add_gene_tags_to_bam.out.CHROM_BAM_BC)

        // cluster_umis(cleanup_headers_3.out.CHROM_BAM_BC_BAI)

        // cleanup_headers_4(cluster_umis.out.bam)

        // combine_chrom_bams(
        //     cleanup_headers_4.outCHROM_BAM_BC_BAI.groupTuple()) // groupBy sample_id?

        // count_cell_gene_umi_reads(combine_chrom_bams.out.BAM_FULLY_TAGGED)

        // umi_gene_saturation(count_cell_gene_umi_reads.out.CELL_UMI_GENE_TSV)

        // construct_expression_matrix(combine_chrom_bams.out.BAM_FULLY_TAGGED)

        // process_expression_matrix(construct_expression_matrix.out.MATRIX_COUNTS_TSV)

        // umap_reduce_expression_matrix(process_expression_matrix.out.MATRIX_PROCESSED_TSV)

        // umap_plot_total_umis(
        //     umap_reduce_expression_matrix.MATRIX_UMAP_TSV
        //     .join(process_expression_matrix.MATRIX_PROCESSED_TSV))

        // TODO: CHANNEL OF GENES
        // genes_to_plot = Channel.fromPath(params.UMAP_PLOT_GENES)
        //     .splitCsv()
        // umap_plot_genes(umap_reduce_expression_matrix.MATRIX_UMAP_TSV
        //     .join(process_expression_matrix.MATRIX_PROCESSED_TSV),
        //     .combine(genes_to_plot))

        // umap_plot_mito_genes(umap_reduce_expression_matrix.MATRIX_UMAP_TSV
        //     .join(process_expression_matrix.MATRIX_PROCESSED_TSV))



}