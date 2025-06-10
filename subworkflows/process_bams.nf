import java.util.ArrayList;

include { merge_and_publish_tsv } from '../modules/local/common'



// Create expression matrices by combining barcode and feature
// tag files. Also outputs the combined tags (per-chrom) to be combined later
process create_matrix {
    label "singlecell"
    cpus 1
    memory "12 GB"
    input:
        tuple val(meta), val(chr), path("features.tsv"), path(read_tags, stageAs: "barcodes.tsv")
    output:
        tuple val(meta), val(chr), path("summary.tsv"), emit: summary
        tuple val(meta), val(chr), path("sa_summary.tsv"), emit: sa_summary
        tuple val(meta), val(chr), val("gene"), path("hdfs/*gene.hdf"), emit: gene
        tuple val(meta), val(chr), val("transcript"), path("hdfs/*transcript.hdf"), emit: transcript
        tuple val(meta), val(chr), path("stats.json"), emit: stats
    script:
    """
    mkdir -p hdfs

    workflow-glue create_matrix \
        ${chr} barcodes.tsv features.tsv \
        --tsv_out summary.tsv \
        --sa_tags_out sa_summary.tsv \
        --hdf_out hdfs \
        --stats stats.json \
        --umi_length ${meta['umi_length']}
    """
}


// Combines multiple expression matrices (e.g. from different chromosomes)
// and calculates summary information on the matrix including UMAPs
process process_matrix {
    label "singlecell"
    cpus  1
    memory "32 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy', pattern: "*{mito,umap,raw,processed}*"
    input:
        tuple val(meta), val(feature), path('inputs/matrix*.hdf')
    output:
        tuple val(meta), val(feature), path("${meta.alias}.${feature}_raw_feature_bc_matrix"), emit: raw
        tuple val(meta), val(feature), path("${meta.alias}.${feature}_processed_feature_bc_matrix"), emit: processed
        tuple val(meta), val(feature), path("${meta.alias}.${feature}_expression_mean_per_cell.tsv"), emit: meancell
        tuple val(meta), val(feature), path("${meta.alias}.${feature}_matrix_stats.tsv"), emit: stats
        // mito per cell makes sense only for feature=gene for now.
        tuple val(meta), val(feature), path("${meta.alias}.gene_expression_mito_per_cell.tsv"), emit: mitocell, optional: true
        tuple val(meta), val(feature), path("${meta.alias}.${feature}_expression_umap*.tsv"), emit: umap
        tuple val(meta), val(feature), path("visium_hd_gene_expression.tsv"), emit: visium_hd, optional: true
    script:
    def mito_prefixes = params.mito_prefix.replaceAll(',', ' ')
    """
    export NUMBA_NUM_THREADS=${task.cpus}
    workflow-glue process_matrix \
        inputs/matrix*.hdf \
        --feature ${feature} \
        --raw "${meta.alias}.${feature}_raw_feature_bc_matrix" \
        --processed "${meta.alias}.${feature}_processed_feature_bc_matrix" \
        --per_cell_mito "${meta.alias}.${feature}_expression_mito_per_cell.tsv" \
        --per_cell_expr "${meta.alias}.${feature}_expression_mean_per_cell.tsv" \
        --umap_tsv "${meta.alias}.${feature}_expression_umap_REPEAT.tsv" \
        --stats "${meta.alias}.${feature}_matrix_stats.tsv" \
        --enable_filtering \
        --min_features $params.matrix_min_genes \
        --min_cells $params.matrix_min_cells \
        --max_mito $params.matrix_max_mito \
        --mito_prefixes $mito_prefixes \
        --norm_count $params.matrix_norm_count \
        --enable_umap \
        --replicates 3
    """
}


// Merge annotated GFFs and transcriptome sequence files
process merge_transcriptome {
    label "singlecell"
    cpus 2
    memory "2 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
            path('fasta/?.fa'),
            path('gffs/?.gff')
    output:
        tuple val(meta),
            path("${meta.alias}.transcriptome.gff.gz"),
            path("${meta.alias}.transcriptome.fa.gz"),
            emit: merged_annotation
    script:
    """
    find fasta/ -name '*.fa' -exec cat {} + \
        | bgzip --threads ${task.cpus} -c  \
        > "${meta.alias}.transcriptome.fa.gz"
    find gffs/ -name '*.gff' -exec cat {} + \
        | grep -v '^#' \
        | bgzip --threads ${task.cpus} -c  \
        > "${meta.alias}.transcriptome.gff.gz"
    """
}


process combine_final_tag_files {
    // Create final per-sample read summaries with information from all stages
    label "singlecell"
    cpus 1
    memory "1 GB"
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              path("tags*.tsv")
    output:
        tuple val(meta),
              path("${meta.alias}.read_summary.tsv")
    script:
    """
    awk 'FNR>1 || NR==1' *.tsv > "${meta.alias}.read_summary.tsv"
    """
}


process umi_gene_saturation {
    label "singlecell"
    cpus 4
    memory {32.GB * task.attempt}
    maxRetries 1
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta),
              path("read_tags.tsv")
    output:
        path("saturation_curves.tsv"),
              emit: saturation_curve
    script:
    """
    export POLARS_MAX_THREADS=$task.cpus

    workflow-glue calc_saturation \
        --output "saturation_curves.tsv" \
        --read_tags read_tags.tsv \
        --sample "${meta.alias}"
    """
}


process pack_images {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta),
              path("images_${meta.alias}/*")
    output:
         tuple val(meta),
              path("images_${meta.alias}")
    script:
    """
    echo packing images
    """
}


process tag_bam {
    label "singlecell"
    cpus 4
    memory {32.GB * task.attempt}
    maxRetries 1
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }

    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    input:
        tuple val(meta),
              path('align.bam'),
              path('align.bam.bai'),
              path('tags/tag_*.tsv'),
              path('sa_tags/sa_tag_*.tsv')
    output:
         tuple val(meta),
               path("${meta.alias}.tagged.bam"),
               path("${meta.alias}.tagged.bam.bai"),
               emit: tagged_bam
    script:
    """
    workflow-glue tag_bam \
        align.bam "${meta.alias}.tagged.bam" tags sa_tags --threads ${task.cpus}
    samtools index "${meta.alias}.tagged.bam"
    """
}


workflow process_bams {
    take:
        merged_bam
        feature_assignments
        sample_annotation
        chr_tags
        read_to_transcript_map  // Rename?
    main:
        // Combine the BAM chunks per-sample
        // Run stringtie per-chrom.
        // Note: this passes in the whole genome BAM but the
        //       .combine() runs this per-chrom such that we get
        //       out reads as fastq per-chrom
        create_matrix(
            feature_assignments
                // Join on [sample meta, chr]
                .join(chr_tags, by: [0, 1]))

        // Aggregate expression matrices to create sparse MEX matrices (https://math.nist.gov/MatrixMarket/formats.html#MMformat)
        // and UMAP TSVs
        process_matrix(
            create_matrix.out.gene.groupTuple(by: [0, 2])
                .map {meta, _chroms, feature, hdfs -> [meta, feature, hdfs.flatten()]}
            .mix(
                create_matrix.out.transcript.groupTuple(by: [0, 2])
                .map {meta, _chroms, feature, hdfs -> [meta, feature, hdfs.flatten()]})
            )

        // TODO: merging the gffs and merging the fasta files is two independent
        //       tasks, they can be done in parallel in two distinct processes.
        merge_transcriptome(
            sample_annotation.groupTuple()
                .join(read_to_transcript_map.groupTuple())
                .map{
                    meta, ann_tr_gff, _chr, tr_fa, _ref_gtf, _str_gff, _fastq ->
                    [meta, tr_fa, ann_tr_gff]})

        // construct per-read summary tables for end user
        // and a tagged bam -- we don't pass final_read_tags here since its
        // advantageous for memory reasons to be able to read the per-chrom
        // tables when iterating over the BAM
        tags_by_sample = create_matrix.out.summary
            .groupTuple()
            .map{meta, _chrs, files -> [meta, files]}
        final_read_tags = combine_final_tag_files(tags_by_sample)

        tag_bam(
            merged_bam.join(tags_by_sample)
            .join(create_matrix.out.sa_summary.groupTuple()
            .map{meta, _chrs, files -> [meta, files]}))
        // UMI saturation curves
        // TODO: this save figures with matplotlib -- just output
        //       data and plot in report with bokeh
        umi_gene_saturation(final_read_tags)

    emit:
        // Emit sperately for use in the report
        // TODO: it shouldn't be the concern of this process what goes in the report
        //       instead just collate everything possible per sample
        final_read_tags = final_read_tags
        tagged_bam = tag_bam.out.tagged_bam
        matrix_stats = process_matrix.out.stats
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        gene_expression = process_matrix.out.processed
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        gene_mean_expression = process_matrix.out.meancell
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        raw_gene_expression = process_matrix.out.raw
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        transcript_mean_expression = process_matrix.out.meancell
            .filter{it[1] == "transcript"}
            .map{it->[it[0], it[2]]}
        mitochondrial_expression = process_matrix.out.mitocell
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
        umap_matrices = process_matrix.out.umap
            .map{it->[it[0], it[2]]}
            .groupTuple(size:2)
            .map{key, files -> [key, files.flatten()]}
        saturation_curves = umi_gene_saturation.out.saturation_curve
            .collectFile(keepHeader: true)
        // per chromosome expression statistics
        expression_stats = create_matrix.out.stats
        // Optional visium HD -
        visium_hd = process_matrix.out.visium_hd
            .filter{it[1] == "gene"}
            .map{it->[it[0], it[2]]}
}
