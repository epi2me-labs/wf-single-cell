include { merge_and_publish_tsv } from '../modules/local/common'
include { cat_tags_by_chrom } from '../modules/local/common'

process generate_whitelist{
    label "singlecell"
    cpus 4
    memory "4 GB"
    publishDir "${params.out_dir}/${meta.alias}",
                mode: 'copy',
                pattern: "*whitelist.tsv"
    input:
        tuple val(meta),
              path("barcodes/?_barcode.tsv")
    output:
        tuple val(meta),
              path("${meta.alias}.whitelist.tsv"),
              emit: whitelist
        path "high_qual_bc_counts.tsv",
              emit: hq_counts
        tuple val(meta),
              path("shortlist_summary.tsv"),
              emit: shortlist_summary
    // TODO: change this to take precomputed, filtered counts from extract_barcodes
    script:
    // It doesn't make sense to do cell count thresholding of the shortlist for visium data.
    // A visium barcode is a tissue coordinate not a cell.
    def no_thresholding_opt = meta.kit.split(':')[0] == 'visium' ? '--no_cell_filter' : ""
    def exp_cells_opt = meta.kit.split(':')[0] != 'visium' ? "--exp_cells ${meta['expected_cells']}" : ""
    def method_opt = params.estimate_cell_count ? "--method quantile" : "--method fixed"
    """
    workflow-glue create_shortlist \
        barcodes "${meta.alias}.whitelist.tsv" shortlist_summary.tsv "${meta.alias}" \
        --counts \
        ${method_opt} \
        ${exp_cells_opt} \
        --counts_out "high_qual_bc_counts.tsv" \
        --threads ${task.cpus} \
        ${no_thresholding_opt}
    """
}


process assign_barcodes{
    label "singlecell"
    cpus 1
    memory "2 GB"
    input:
         tuple val(meta),
               path("whitelist.tsv"),
               path("extract_barcodes.tsv")
    output:
        tuple val(meta),
              path("bc_assign_counts.tsv"),
              emit: chrom_assigned_barcode_counts
        tuple val(meta),
              path("extract_barcodes_with_bc.tsv"),
              emit: tags
        tuple val(meta),
              path("summary.tsv"),
              emit: summary
    script:
    """
    workflow-glue assign_barcodes \
        whitelist.tsv extract_barcodes.tsv \
        extract_barcodes_with_bc.tsv bc_assign_counts.tsv summary.tsv \
        --max_ed ${params.barcode_max_ed} \
        --min_ed_diff ${params.barcode_min_ed_diff} \
        --use_kmer_index
    """
}



workflow correct_10x_barcodes {
    take: 
        extracted_barcodes
        high_qual_bc_counts
    main:


        generate_whitelist(high_qual_bc_counts)

        // TODO: this process really has no business being here. It should be
        //       moved into main.nf as an aggregation across all the chunks
        //       in extracted_barcodes. It takes a long time per-chunk so should
        //       be left as parallel across chunks.
        assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extracted_barcodes)
            .map {it ->
                def meta = it[0][0]
                def whitelist = it[0][1]
                def barcodes = it[1][1]
                [meta, whitelist, barcodes]})

        merge_and_publish_tsv(
            assign_barcodes.out.summary
                .concat(generate_whitelist.out.shortlist_summary)
                .groupTuple(),
            'bc_assignment_summary.tsv')

        // Combine the tag chunks to per chrom chunks and emit [meta, chr, tags]
        chr_tags = cat_tags_by_chrom(assign_barcodes.out.tags.groupTuple())
            .transpose()
            .map {meta, file -> [meta, file.baseName, file]}

    emit:
        chr_tags = chr_tags
        white_list = generate_whitelist.out.whitelist
        hq_bc_counts = generate_whitelist.out.hq_counts.collectFile(keepHeader: true)
}