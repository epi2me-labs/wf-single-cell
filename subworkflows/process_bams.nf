import java.util.ArrayList;

include { merge_and_publish_tsv } from '../modules/local/common'

process split_gtf_by_chroms {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path("ref.gtf")
    output:
        path("*"), emit: chrom_gtf
    script:
    """
    gawk '/^[^#]/ {print>\$1".gtf"}' ref.gtf
    """
}


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


process merge_bams {
    // Combine all BAMs derived from the initial chunking into per sample files
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
            path('bams/*aln.bam'),
            path('bams/*aln.bam.bai')
    output:
        tuple val(meta),
              path("merged.sorted.bam"),
              path("merged.sorted.bam.bai"),
              emit: merged_bam
    script:
    """
    samtools merge -@ ${task.cpus -1} --write-index -o "merged.sorted.bam##idx##merged.sorted.bam.bai" bams/*.bam
    """
}


process cat_tags_by_chrom {
    // Merge per-chunk tags to create per-chromosome tags
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path('tags/*tags.tsv')
    output:
         tuple val(meta),
              path("chr_tags/*"),
              emit: merged_tags

    script:
    """
    mkdir chr_tags
    # Find the chr column number
    files=(tags/*)
    chr_col=\$(awk -v RS='\t' '/chr/{print NR; exit}' "\${files[0]}")

    # merge the tags TSVs, keep header from first file and split entries by chromosome
    awk -F'\t' -v chr_col=\$chr_col 'FNR==1{hdr=\$0; next} \
    {if (!seen[\$chr_col]++) \
        print hdr>"chr_tags/"\$chr_col".tsv"; \
        print>"chr_tags/"\$chr_col".tsv"}' tags/*

    # Sort by corrected barcode to allow chunked reading later.
    for file in chr_tags/*.tsv;
    do
        csvtk sort -t --keys CB "\$file" -o tmp.tsv;
        mv tmp.tsv "\$file";
    done
    """
}


process stringtie {
    label "singlecell"
    cpus params.threads
    // Memory usage for this process is usually less than 3GB, but some cases it may go over this.
    memory { 3.GB * task.attempt }
    maxRetries 3
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        path 'ref_genome.fa'
        path 'ref_genome.fa.fai'
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai"),
              val(chr),
              path("chr.gtf")

    output:
        tuple val(meta),
              val(chr),
              path("transcriptome.fa"),
              path("chr.gtf"),
              path("stringtie.gff"),
              path("reads.fastq.gz"),
              emit: read_tr_map
    script:
    """
    # Add chromosome label (-l) to generated transcripts
    # so we don't get name collisions during file merge later
    samtools view -h align.bam ${chr}  \
        | tee >(
            stringtie -L ${params.stringtie_opts} -p ${task.cpus} \
                -G chr.gtf -l "${chr}.stringtie" -o "stringtie.gff" - ) \
        | samtools fastq \
        | bgzip --threads 2 -c > reads.fastq.gz
    # Get transcriptome sequence
    gffread -g ref_genome.fa -w "transcriptome.fa" "stringtie.gff"
    """
}


process align_to_transcriptome {
    label "singlecell"
    cpus params.threads
    memory "31 GB"
    input:
        tuple val(meta),
              val(chr),
              path('transcriptome.fa'),
              path('chr.gtf'),
              path('stringtie.gff'),
              path("reads.fq.gz")
    output:
        tuple val(meta),
              val(chr),
              path("chr.gtf"),
              path("tr_align.bam"),
              path('stringtie.gff'),
              emit: read_tr_map
    script:
    def view_threads = 1
    def sort_threads = 3
    def mm2_threads = Math.max(task.cpus - view_threads - sort_threads, 4)
    """
    minimap2 -ax map-ont \
        --cap-kalloc 100m --cap-sw-mem 50m \
        --end-bonus 10 -p 0.9 -N 3 -t $mm2_threads \
        transcriptome.fa reads.fq.gz \
    | samtools view -h -@ $view_threads -b -F 2052 - \
    | samtools sort -n -@ $sort_threads --no-PG - > tr_align.bam
    """
}


process assign_features {
    label "singlecell"
    cpus 1
    // This step is performed per-chromosome. The tags file per chrom can vary
    // quite widely in size. We don't have a fixed memory size here in order
    // to get better parallelism on single-host setups.
    memory { 1.0.GB.toBytes() + (tags.size() * 2 ) }
    input:
        tuple val(meta),
              val(chr),
              path("chr.gtf"),
              path("tr_align.bam"),
              path("stringtie.gff"),
              path(tags, stageAs: "tags.tsv")
    output:
        tuple val(meta),
              val(chr),
              path("feature_assigns.tsv"),
              emit: feature_assigns
        tuple val(meta),
              path("gffcompare.annotated.gtf"),
              emit: annotation
    script:
    """
    # gffcomapre maps transcript reference IDs to query transcripts.
    gffcompare -o gffcompare -r chr.gtf stringtie.gff

    workflow-glue assign_features \
        tr_align.bam \
        gffcompare.stringtie.gff.tmap \
        chr.gtf \
        tags.tsv \
        feature_assigns.tsv \
        --min_mapq ${params.gene_assigns_minqv}
    """
}


// Create expression matrices by combining barcode and feature
// tag files. Also outputs the combined tags (per-chrom) to be combined later
process create_matrix {
    label "singlecell"
    cpus 1
    memory "7 GB"
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
    memory "16 GB"
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
    memory "31 GB"
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
    memory "31 GB"
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
        bam
        extracted_barcodes
        high_qual_bc_counts
        gtf
        ref_genome_fasta
        ref_genome_idx
    main:
        // Split the GTF by chromosome
        chr_gtf = split_gtf_by_chroms(gtf)
            .flatten()
            .map {fname -> tuple(fname.baseName, fname)} // [chr, gtf]

        generate_whitelist(high_qual_bc_counts)

        // TODO: this process really has no business being here. It should be
        //       moved into main.nf as an aggregation across all the chunks
        //       in extracted_barcodes. It takes a long time per-chunk so should
        //       be left as parallel across chunks.
        assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extracted_barcodes)
            .map {it ->
                meta = it[0][0]
                whitelist = it[0][1]
                barcodes = it[1][1]
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

        // Combine the BAM chunks per-sample
        merge_bams(bam.groupTuple())

        // Run stringtie per-chrom.
        // Note: this passes in the whole genome BAM but the
        //       .combine() runs this per-chrom such that we get
        //       out reads as fastq per-chrom
        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            merge_bams.out.merged_bam
                .combine(chr_gtf))

        // TODO: We're likely to change this to use bambu and avoid using
        //       stringtie altogether. However note that the next three steps
        //       are a strict linear pipeline and should be combined into one
        //       process to avoid staging of files between processes. Note further
        //       that it would be trivial to combine the assign_features and
        //       and create_matrix steps into a single program to avoid writing
        //       any intermediate files whatsoever.
        align_to_transcriptome(stringtie.out.read_tr_map)

        assign_features(
            align_to_transcriptome.out.read_tr_map
                .join(chr_tags, by: [0, 1]))

        create_matrix(
            assign_features.out.feature_assigns
                // Join on [sample meta, chr]
                .join(chr_tags, by: [0, 1]))

        // Aggregate expression matrices to create sparse MEX matrices (https://math.nist.gov/MatrixMarket/formats.html#MMformat)
        // and UMAP TSVs
        process_matrix(
            create_matrix.out.gene.groupTuple(by: [0, 2])
                .map {meta, chroms, feature, hdfs -> [meta, feature, hdfs.flatten()]}
            .mix(
                create_matrix.out.transcript.groupTuple(by: [0, 2])
                .map {meta, chroms, feature, hdfs -> [meta, feature, hdfs.flatten()]})
            )

        // TODO: merging the gffs and merging the fasta files is two independent
        //       tasks, they can be done in parallel in two distinct processes.
        merge_transcriptome(
            assign_features.out.annotation.groupTuple()
                .join(stringtie.out.read_tr_map.groupTuple())
                .map{
                    meta, ann_tr_gff, chr, tr_fa, ref_gtf, str_gff, fastq ->
                    [meta, tr_fa, ann_tr_gff]})

        // construct per-read summary tables for end user
        // and a tagged bam -- we don't pass final_read_tags here since its
        // advantageous for memory reasons to be able to read the per-chrom
        // tables when iterating over the BAM
        tags_by_sample = create_matrix.out.summary
            .groupTuple()
            .map{meta, chrs, files -> [meta, files]}
        final_read_tags = combine_final_tag_files(tags_by_sample)

        tag_bam(
            merge_bams.out.join(tags_by_sample)
            .join(create_matrix.out.sa_summary.groupTuple()
            .map{meta, chrs, files -> [meta, files]}))
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
        white_list = generate_whitelist.out.whitelist
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
        hq_bc_counts = generate_whitelist.out.hq_counts.collectFile(keepHeader: true)
        // per chromosome expression statistics
        expression_stats = create_matrix.out.stats
}
