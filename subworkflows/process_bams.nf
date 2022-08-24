import java.util.ArrayList;

process get_kit_info {
    label "singlecell"
    
    input:
        path kit_config
        path sc_sample_sheet
              
    output:
        path 'sample_kit_info.csv', emit: kit_info
    
    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    sample_df = pd.read_csv("${sc_sample_sheet}", sep=",", comment="#").set_index("sample_id", drop=True)

    kit_df = pd.read_csv("$kit_config", sep=",", comment="#")


    records = []
    for sample_id, row in sample_df.iterrows():
        kit_name = sample_df.loc[sample_id, "kit_name"]
        kit_version = sample_df.loc[sample_id, "kit_version"]

        # Get barcode length based on the kit_name and kit_version specified for this sample_id.
        rows = (kit_df["kit_name"] == kit_name) & (kit_df["kit_version"] == kit_version)
        barcode_length = kit_df.loc[rows, "barcode_length"].values[0]

        # Get the appropriate cell barcode longlist based on the kit_name specified for this sample_id.
        if kit_name == "3prime":
            long_list = "3M-february-2018.txt.gz"
        elif kit_name == "5prime":
            long_list = "737K-august-2016.txt.gz"
        elif kit_name == "multiome":
            long_list = "737K-arc-v1.txt.gz"
        else:
            raise Exception("Encountered an unexpected kit_name in samples.csv")

        # Get UMI length based on the kit_name and kit_version specified for this sample_id.
        umi_length = kit_df.loc[rows, "umi_length"].values[0]
        records.append([sample_id, kit_name, kit_version, barcode_length, umi_length, long_list])
    df_out = pd.DataFrame.from_records(records)
    df_out.columns = ["sample_id", "kit_name", "kit_version", "barcode_length", "umi_length", "long_list"]
    df_out.to_csv('sample_kit_info.csv', index=False)
    """
}


process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    
    cpus params.max_threads
    input:
        tuple val(sample_id), 
              val(kit_name),
              val(kit_version),
              val(barcode_length),
              val(umi_length),
              val(bc_long_list),
              path(bam_sort),
              path(bam_sort_idx)
        path bc_longlist_dir

    output:
        tuple val(sample_id), path("*.bam"), emit: bam_bc_uncorr_tmp
        tuple val(sample_id), path("*.tsv"), emit: barcode_counts
    """
    extract_barcode.py \
    $bam_sort ${bc_longlist_dir}/${bc_long_list}\
    -t $task.cpus \
    --kit $kit_name \
    --adapter1_suff_length $params.barcode_adapter1_suff_length \
    --min_barcode_qv $params.barcode_min_quality \
    --barcode_length $barcode_length \
    --umi_length $umi_length \
    --output_bam "tmp.bam" \
    --output_barcodes "${sample_id}.uncorrected_bc_counts.tsv";
    """
}

process cleanup_headers_1 {
    label "singlecell"

    cpus params.max_threads
    
    input:
         tuple val(sample_id), path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bam.bai"), 
        emit: bam_bc_uncorr
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam > "${sample_id}.bc_extract.sorted.bam"
    samtools index "${sample_id}.bc_extract.sorted.bam"
    """
}

process generate_whitelist{
    label "singlecell"
    
    cpus params.max_threads
    input:
        tuple val(sample_id), 
              path(counts),
              val(expected_cells)
    output:
        tuple val(sample_id), path("*whitelist.tsv"), emit: whitelist
        tuple val(sample_id), path("*kneeplot.png"), emit: kneeplot
     def kneeflags = params.barcode_kneeplot_flags ?  params.barcode_kneeplot_flags : ''
    
    """
    knee_plot.py \
        ${kneeflags} \
        --exp_cells $expected_cells \
        --output_whitelist "${sample_id}.whitelist.tsv" \
        --output_plot "${sample_id}.kneeplot.png" $counts
    """
}

process split_bam_by_chroms{
    label "singlecell"
    
    cpus params.max_threads
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        tuple val(sample_id), 
              path("splits/*.bam"),
              path("splits/*.bai"), emit: bam
    """
    split_bam_by_chroms.py -t ${task.cpus} --output_dir splits $bam
    """
}


process assign_barcodes{
    label "singlecell"
    
    input:
         tuple val(sample_id), 
               val(kit_name),
               val(kit_version),
               val(barcode_length),
               val(umi_length),
               val(bc_long_list),
               path(whitelist),
               val(_), // Redundant sample_id: remove
               val(chr), 
               path(bam),
               path(bai)
    output:
        tuple val(sample_id), 
            val(chr),
            path("tmp.bam"),
            emit: chrom_bam_bc
        tuple val(sample_id),
              val(chr),
              path("*.bc_assign_counts.tsv"), 
              emit: chrom_assigned_barcode_counts
    """
    assign_barcodes.py -t 1 \
        --output_bam tmp.bam \
        --output_counts ${sample_id}_${chr}.bc_assign_counts.tsv \
        --max_ed $params.barcode_max_ed \
        --min_ed_diff $params.barcode_min_ed_diff \
        --kit $kit_name \
        --adapter1_suff_length $params.barcode_adapter1_suff_length \
        --barcode_length $barcode_length \
        --umi_length $umi_length \
        $bam $whitelist
"""
}

process cleanup_headers_2 {
    label "singlecell"

    cpus params.max_threads
    
    input:
         tuple val(sample_id), val(chr), path(bam)
    output:
        tuple val(sample_id), 
              val(chr), 
              path("*.bam"), 
              path("*.bam.bai"), 
        emit: chrom_bam_bc_bai
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        > "${sample_id}_${chr}.bc_assign.bam";
    samtools index "${sample_id}_${chr}.bc_assign.bam"
    """
}

process bam_to_bed {
    label "singlecell"
    
    input:
        tuple val(chr), //emit chr first for doing cross on gtfs 
              val(sample_id),
              path(bam),
              path(bai) 
    output:
        tuple val(sample_id), 
              val(chr),
              path('*bc_assign.bed'), 
              emit: chrom_bed_bc
    """
    bedtools bamtobed -i $bam > "${sample_id}_bc_assign.bed"
    """
}

process split_gtf_by_chroms {
    label "singlecell"
    input:
        path(gtf)
    output:
        path("*"), emit: chrom_gtf
    """
    awk '/^[^#]/ {print>\$1}' $gtf 
    """
}   

process assign_genes {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              val(chr),
              path(chrom_bed_bc),
              path(chrom_gtf)
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_${chr}.read.gene_assigns.tsv"),
              emit: chrom_tsv_gene_assigns
    """
    assign_genes.py \
    --output ${sample_id}_${chr}.read.gene_assigns.tsv \
    $chrom_bed_bc $chrom_gtf
    """
}

process add_gene_tags_to_bam {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              val(chr),
              path(chrom_bam_bc),
              path(chrom_bam_bc_bai),
              path(chrom_tsv_gene_assigns)
    output:
        tuple val(sample_id),
              val(chr),
              path("tmp.bam"), 
              emit: chrom_bam_bc_gene_tmp
    """
    add_gene_tags.py \
        --output tmp.bam \
        $chrom_bam_bc $chrom_tsv_gene_assigns
    """
}

process cleanup_headers_3 {
    label "singlecell"
    cpus params.max_threads
    
    input:
         tuple val(sample_id), 
               val(chr), 
               path(bam)
    output:
        tuple val(sample_id), 
              val(chr),
              path("*.bam"), 
              path("*.bam.bai"), 
              emit: chrom_bam_bc_bai
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        > ${sample_id}_${chr}_bc_assign.gene.bam;
    samtools index ${sample_id}_${chr}_bc_assign.gene.bam
    """
}

process cluster_umis {
    label "singlecell"
    
    cpus params.umi_cluster_max_threads
    input:
        tuple val(sample_id),
              val(chr),
              path(bam),
              path(bai)
    output:
         tuple val(sample_id),
              val(chr),
              path("tmp.bam"),
              emit: bam
    """
    cluster_umis.py \
    $bam \
    --threads $task.cpus \
    --ref_interval $params.umi_genomic_interval \
    --cell_gene_max_reads $params.umi_cell_gene_max_reads \
    --output tmp.bam 
    """
}

process cleanup_headers_4{
    label "singlecell"
    
    cpus params.max_threads
    input:
         tuple val(sample_id), 
               val(chr), 
               path(bam)
    output:
        tuple val(sample_id), path("*.bam"), path("*.bam.bai"), 
        emit: chrom_bam_bc_bai
    """
    samtools reheader --no-PG -c 'grep -v ^@PG' $bam \
        >  ${sample_id}_${chr}.tagged.bam;
    samtools index ${sample_id}_${chr}.tagged.bam
    """
}

process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "singlecell"
    
    input:
        tuple val(sample_id), 
              path(bams),
              path(bais)
    output:
        tuple val(sample_id), 
              path("*tagged.sorted.bam"), 
              path("*tagged.sorted.bam.bai"),
              emit: bam_fully_tagged
    """
    samtools merge -o "${sample_id}.tagged.sorted.bam" $bams; 
    samtools index "${sample_id}.tagged.sorted.bam";
    """
}

process count_cell_gene_umi_reads {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(bam),
              path(bai)
    output:
        tuple val(sample_id),
              path("${sample_id}_cell_umi_gene.tsv"),
              emit: cell_umi_gene_tsv
    """
    cell_umi_gene_table.py \
        --output  ${sample_id}_cell_umi_gene.tsv $bam
    """
}

process umi_gene_saturation {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(cell_umi_gene_tsv)
    output:
        tuple val(sample_id),
              path("*saturation_curves.png")
    """
    calc_saturation.py \
        --output ${sample_id}.saturation_curves.png \
        $cell_umi_gene_tsv
    """
}

process construct_expression_matrix {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(bam),
              path(bai)
    output:
        tuple val(sample_id), 
              path("*gene_expression.counts.tsv"), 
              emit: matrix_counts_tsv
    """
    gene_expression.py \
        --output ${sample_id}.gene_expression.counts.tsv $bam
    """
}

process process_expression_matrix {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(matrix_counts_tsv)
    output:
        tuple val(sample_id), 
              path("*gene_expression.processed.tsv"),
              emit: matrix_processed_tsv
        tuple val(sample_id),
              path("*gene_expression.mito.tsv"),
              emit: matrix_mito_tsv
    """
    process_matrix.py \
    --min_genes $params.matrix_min_genes \
    --min_cells $params.matrix_min_cells \
    --max_mito $params.matrix_max_mito \
    --norm_count $params.matrix_norm_count \
    --output ${sample_id}.gene_expression.processed.tsv \
    --mito_output ${sample_id}.gene_expression.mito.tsv \
    $matrix_counts_tsv
    """
}

process umap_reduce_expression_matrix {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(matrix_processed_tsv)
    output:
         tuple val(sample_id),
              path("*gene_expression.umap.tsv"), 
              emit: matrix_umap_tsv
    """
    umap_reduce.py \
        --output ${sample_id}.gene_expression.umap.tsv \
        $matrix_processed_tsv
    """
}


process umap_plot_total_umis {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_processed_tsv)
    output:
          tuple val(sample_id),
              path("*umap.total.png"), 
              emit: matrix_umap_plot_total
    """
    plot_umap.py \
        --output ${sample_id}.umap.total.png \
        $matrix_umap_tsv $matrix_processed_tsv
    """
}

process umap_plot_genes {
    // TODO: make a channle of input genes for thes process
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_processed_tsv),
              val(gene)
    output:
        tuple val(sample_id),
              path("*umap.gene.${gene}.png"), 
              emit: matrix_umap_plot_gene
    """
    plot_umap.py \
        --gene $gene \
        --output ${sample_id}.umap.gene.${gene}.png \
        $matrix_umap_tsv $matrix_processed_tsv
    """
}

process umap_plot_mito_genes {
    label "singlecell"
    
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_mito_tsv)
    output:
        tuple val(sample_id),
              path("*umap.mitochondrial.png"), 
              emit: matrix_umap_plot_mito
    """
    plot_umap.py \
        --mito_genes \
        --output ${sample_id}.umap.mitochondrial.png \
        $matrix_umap_tsv \
        $matrix_mito_tsv
    """
}
    
workflow process_bams {
    take:
        bam
        bam_idx
        sc_sample_sheet
        kit_config
        gtf
        umap_genes
        bc_longlist_dir
        sample_kits
   
    main:

        get_kit_info(
            kit_config,
            sc_sample_sheet)

        extract_barcodes(
            get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(bam).join(bam_idx),
            bc_longlist_dir)
        
        cleanup_headers_1(
            extract_barcodes.out.bam_bc_uncorr_tmp
        )
        
        split_bam_by_chroms(
            cleanup_headers_1.out.bam_bc_uncorr
        )

        generate_whitelist(
            extract_barcodes.out.barcode_counts
            .join(sample_kits)
            .map{it -> tuple(it[0], it[1], it[4] )}
        )
        
        // Extract chr from filename and add to tuple to give: 
        // [sample_id, chr, bam, bai]
        bam_bai_chromes = split_bam_by_chroms.out.bam.map({it ->
            sbi = []
            if (it[1].getClass() == java.util.ArrayList){
                // Multiple chroms:
                // [sample_id, [bam1, bam2], [bai1, bai2]]
                for (i=0; i<it[1].size(); i++) {
                    chr = it[1][i].toString().tokenize('/')[-1] - '.sorted.bam'
                    sbi.add(tuple(it[0], chr, it[1][i], it[2][i]))
                }
            }
            else{
                println(it[1].getClass())
                // Only single chrom so we have:
                // [sample_id, bam, bai]
                chr = it[1].toString().tokenize('/')[-1] - '.sorted.bam'
                sbi.add(tuple(it[0], chr, it[1], it[2])) 
            }
            return sbi
            }).flatMap(it-> it)
        
        // merge kit info to bams
        chr_bam_kit = get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(generate_whitelist.out.whitelist)
            .cross(bam_bai_chromes).map({it -> it.flatten()})

        assign_barcodes(chr_bam_kit)

        cleanup_headers_2(assign_barcodes.out.chrom_bam_bc)

        bam_to_bed(cleanup_headers_2.out.chrom_bam_bc_bai)

        chr_gtf = split_gtf_by_chroms(gtf)
        .flatten()
        .map {file -> 
            // create [chr, gtf]
             tuple(file.toString().tokenize('/')[-1], file)}
        
        // combine all chr bams with chr gtfs
        chr_bams_gtf = chr_gtf.cross(
            bam_to_bed.out.chrom_bed_bc)
            .map({it ->
            //  rejig the tuple to [sample_id, chr, bed, gtf]
             tuple(it[1][1], it[0][0], it[1][2], it[0][1])})

         assign_genes(chr_bams_gtf)

         add_gene_tags_to_bam(
             cleanup_headers_2.out.chrom_bam_bc_bai
              // join on sample_id + chr
             .join(assign_genes.out.chrom_tsv_gene_assigns, by:[0, 1]))

         cleanup_headers_3(add_gene_tags_to_bam.out.chrom_bam_bc_gene_tmp)

         cluster_umis(cleanup_headers_3.out.chrom_bam_bc_bai)

         cleanup_headers_4(cluster_umis.out.bam)

         // group by sample_id
         combine_chrom_bams(
             cleanup_headers_4.out.chrom_bam_bc_bai.groupTuple())

         count_cell_gene_umi_reads(combine_chrom_bams.out.bam_fully_tagged)

         umi_gene_saturation(count_cell_gene_umi_reads.out.cell_umi_gene_tsv)

         construct_expression_matrix(combine_chrom_bams.out.bam_fully_tagged)

         process_expression_matrix(construct_expression_matrix.out.matrix_counts_tsv)

         umap_reduce_expression_matrix(process_expression_matrix.out.matrix_processed_tsv)

         umap_plot_total_umis(
             umap_reduce_expression_matrix.out.matrix_umap_tsv
             .join(process_expression_matrix.out.matrix_processed_tsv))

         genes_to_plot = Channel.fromPath(umap_genes)
             .splitCsv()
        
         umap_plot_genes(
             umap_reduce_expression_matrix.out.matrix_umap_tsv
             .join(process_expression_matrix.out.matrix_processed_tsv)
             .combine(genes_to_plot))

         umap_plot_mito_genes(
            umap_reduce_expression_matrix.out.matrix_umap_tsv
            .join(process_expression_matrix.out.matrix_mito_tsv))
        
     emit:
         umap_plots = umap_plot_genes.out.groupTuple()
            .flatMap({it -> 
                l = []
                for (i=0; i<it[1].size(); i++)
                    l.add(tuple(it[0], it[1][i]))
                return l
                })
                .concat(umap_plot_total_umis.out)
                .concat(umap_plot_mito_genes.out)
                .groupTuple()
                
        results = umi_gene_saturation.out
             .join(construct_expression_matrix.out)
             .join(process_expression_matrix.out.matrix_processed_tsv)
             .join(process_expression_matrix.out.matrix_mito_tsv)
             .join(cleanup_headers_1.out.bam_bc_uncorr)
             .join(generate_whitelist.out.whitelist)
             .join(generate_whitelist.out.kneeplot)
             .join(combine_chrom_bams.out.bam_fully_tagged)
             .join(extract_barcodes.out.barcode_counts)
             .join(umap_reduce_expression_matrix.out.matrix_umap_tsv)
}