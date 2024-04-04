import java.util.ArrayList;

include {
    construct_expression_matrix as construct_gene_expression;
    construct_expression_matrix as construct_transcript_expression;
} from '../modules/local/common'


process get_contigs {
    label "singlecell"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta),
              path('sample.bam'),
              path('sample.bam.bai')

    output:
        tuple path("${meta.alias}_contigs"),
            val(meta),
            emit: contigs
    """
    samtools idxstats sample.bam \
        | gawk '/^[^*]/{print\$1}' \
        | gawk NF > "${meta.alias}_contigs"
    """
}


process generate_whitelist{
    label "singlecell"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta),
              path("barcodes/?_barcode.tsv")
        path("barcode_longlist_dir")
    output:
        tuple val(meta),
              path("*whitelist.tsv"), 
              emit: whitelist
        tuple val(meta),
              path("*kneeplot.png"), 
              emit: kneeplot
        tuple val(meta),
              path("${meta.alias}.uncorrected_bc_counts.tsv"), 
              emit: uncorrected_bc_counts
    """
    workflow-glue knee_plot \
        --barcodes_dir barcodes/ \
        --long_list "barcode_longlist_dir/${meta['bc_long_list']}" \
        --exp_cells ${meta['expected_cells']} \
        --output_whitelist "${meta.alias}.whitelist.tsv" \
        --output_plot "${meta.alias}.kneeplot.png" \
        --output_uncorrected_barcodes "${meta.alias}.uncorrected_bc_counts.tsv"
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
    """
    workflow-glue assign_barcodes \
        --output_tags extract_barcodes_with_bc.tsv \
        --output_counts bc_assign_counts.tsv \
        --max_ed $params.barcode_max_ed \
        --min_ed_diff $params.barcode_min_ed_diff \
        --extract_barcode_tags extract_barcodes.tsv \
        --whitelist whitelist.tsv
    """
}


process split_gtf_by_chroms {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path("ref.gtf")
    output:
        path("*"), emit: chrom_gtf
    """
    gawk '/^[^#]/ {print>\$1".gtf"}' ref.gtf 
    """
}   


process cluster_umis {
    label "singlecell"
    cpus 1
    // Benchmarking showed that memory usage was ~ 15x the size of read_tags input.
    // Set a minimum memory requirement of 1.0GB to allow for overhead.
    memory {1.0.GB.toBytes()  + (read_tags.size() * 20) }
    input:
        tuple val(meta),
              val(chr),
              path("chrom_feature_assigns.tsv"),
              path(read_tags, stageAs: "read_tags.tsv")
    output:
        tuple val(meta),
              val(chr),
              path("${meta.alias}_${chr}.read_tags.tsv"),
              emit: read_tags  // For BAM tagging
        tuple val(meta),
              path("${meta.alias}_${chr}.final_tags.tsv"),
              emit: final_read_tags  // For user output
    """
    workflow-glue cluster_umis \
        --chrom ${chr} \
        --feature_assigns chrom_feature_assigns.tsv \
        --read_tags read_tags.tsv \
        --output_read_tags "${meta.alias}_${chr}.read_tags.tsv" \
        --workflow_output "${meta.alias}_${chr}.final_tags.tsv"
    """
}


process tag_bams {
    label "singlecell"
    cpus 4
    memory "16 GB"
    input:
        tuple val(meta),
              path("align.bam"),
              path("align.bam.bai"),
              val(chr),
              path('tags.tsv')
    output:
         tuple val(meta),
              path("${meta.alias}.${chr}.tagged.bam"),
              path("${meta.alias}.${chr}.tagged.bam.bai"),
              emit: tagged_bam
    script:
    """
    workflow-glue tag_bam \
        --in_bam align.bam \
        --tags tags.tsv \
        --out_bam ${meta.alias}.${chr}.tagged.bam \
        --chrom ${chr}

    samtools index -@ ${task.cpus} ${meta.alias}.${chr}.tagged.bam
    """
}


process combine_tag_files {
    // A file named 'read_tags' is a required output.
    // collectFile's 'name' argument does not work when being applied to
    // a channel that returns tuples. It groups and names according to the
    // first element of the tuple. Hence this process.
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta),
              path("tags*.tsv")
    output:
        tuple val(meta),
              path("${meta.alias}_read_tags.tsv")
    """
    awk 'FNR>1 || NR==1' *.tsv > "${meta.alias}_read_tags.tsv"
    """
}


process combine_final_tag_files {
    // Combine the final
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta),
              path("tags*.tsv")
    output:
        tuple val(meta),
              path("${meta.alias}.read_tags.tsv")
    """
    awk 'FNR>1 || NR==1' *.tsv > "${meta.alias}.read_tags.tsv"
    """
}


process combine_bams_and_tags {
    // Merge all BAM and tags files chunks
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path('bams/*aln.bam'),
              path('bams/*.aln.bam.bai'),
              path('tags/*tags.tsv')
    output:
        tuple val(meta),
              path("*tagged.sorted.bam"), 
              path("*tagged.sorted.bam.bai"),
              emit: merged_bam
        tuple val(meta),
              path("chr_tags/*"),
              emit: merged_tags
    """
    samtools cat -b <(find bams -name '*aln.bam') \
    | samtools sort - -@ ${task.cpus} --write-index \
        -o "${meta.alias}.tagged.sorted.bam##idx##${meta.alias}.tagged.sorted.bam.bai"

    mkdir chr_tags
    # Find the chr column number
    files=(tags/*)
    chr_col=\$(awk -v RS='\t' '/chr/{print NR; exit}' "\${files[0]}")

    # merge the tags TSVs, keep header from first file and split entries by chromosome
    awk -F'\t' -v chr_col=\$chr_col 'FNR==1{hdr=\$0; next} \
    {if (!seen[\$chr_col]++) \
        print hdr>"chr_tags/"\$chr_col".tsv"; \
        print>"chr_tags/"\$chr_col".tsv"}' tags/*
    """
}

process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path(chrom_bams),
              path('chrom.bam.bai')
    output:
        tuple val(meta),
              path("*tagged.sorted.bam"),
              path("*tagged.sorted.bam.bai"),
              emit: bam_fully_tagged
    """
    samtools merge -@ ${task.cpus - 1} --write-index \
        -o "${meta.alias}.tagged.sorted.bam##idx##${meta.alias}.tagged.sorted.bam.bai" ${chrom_bams};
    """
}


process stringtie {
    label "singlecell"
    cpus params.threads
    // Memory usage for this process is usually less than 3GB, but some cases it may go over this.
    memory = { 3.GB * task.attempt }
    maxRetries = 3
    errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
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
              path("${meta.alias}.transcriptome.fa"),
              path("chr.gtf"),
              path("${meta.alias}.stringtie.gff"),
              path("reads.fastq"),
              emit: read_tr_map
    script:
    """
    # Add chromosome label (-l) to generated transcripts
    # so we don't get name collisions during file merge later
    samtools view -h align.bam ${chr}  \
         | tee >(stringtie -L ${params.stringtie_opts} -p ${task.cpus} -G chr.gtf -l "${chr}.stringtie" \
             -o "${meta.alias}.stringtie.gff" - ) \
         | samtools fastq > reads.fastq
    # Get transcriptome sequence
    gffread -g ref_genome.fa -w "${meta.alias}.transcriptome.fa" "${meta.alias}.stringtie.gff"
    """
}


process align_to_transcriptome {
    label "singlecell"
    cpus params.threads
    memory = "32 GB"
    input:
        tuple val(meta),
              val(chr),
              path('transcriptome.fa'),
              path('chr.gtf'),
              path('stringtie.gff'),
              path("reads.fq")
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
        transcriptome.fa reads.fq \
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
              path('stringtie.gff'),
              path(tags, stageAs: 'tags.tsv')
    output:
        tuple val(meta),
              val(chr),
              path("${meta.alias}.${chr}.feature_assigns.tsv"),
              emit: feature_assigns
        tuple val(meta),
              path("gffcompare.annotated.gtf"),
              emit: annotation
    """
    # gffcomapre maps transcript reference IDs to query transcripts.
    gffcompare -o gffcompare -r chr.gtf stringtie.gff

    workflow-glue assign_features \
        --transcriptome_bam tr_align.bam \
        --gffcompare_tmap gffcompare.stringtie.gff.tmap \
        --gtf chr.gtf \
        --tags tags.tsv \
        --output "${meta.alias}.${chr}.feature_assigns.tsv" \
        --min_mapq ${params.gene_assigns_minqv}
    """
}

process umi_gene_saturation {
    label "singlecell"
    cpus 4
    memory "32 GB"
    input:
        tuple val(meta),
              path("read_tags.tsv")
    output:
        tuple val(meta),
              path("*saturation_curves.png"),
              emit: saturation_curve
    """
    export POLARS_MAX_THREADS=$task.cpus

    workflow-glue calc_saturation \
        --output "${meta.alias}.saturation_curves.png" \
        --read_tags read_tags.tsv
    """
}


process umap_reduce_expression_matrix {
    label "singlecell"
    cpus 2
    // Most runs will use less than 10GB memory, but large numbers of cells (above 15K)
    // may lead to memory usage over 20GB. Max here is 32GB 
    memory { 8.GB * task.attempt }
    maxRetries 3
    errorStrategy {task.exitStatus in [137,140] ? 'retry' : 'terminate'}
    input:
        tuple val(repeat_num),
              val(meta),
              val(data_type),
              path(matrix)
    output:
        tuple val(meta),
                path("${data_type}_umap_${repeat_num}.tsv"),
                emit: matrix_umap_tsv
    """
    workflow-glue umap_reduce \
        --output ${data_type}_umap_${repeat_num}.tsv \
        ${matrix}
    """
}


process merge_transcriptome {
    // Merge the annotated GFFs and transcriptome sequence files
    label "singlecell"
    cpus 1
    memory "2GB"
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
    # Concatenate transcriptome files, remove comments (from gff) and compress
    find fasta/ -name '*.fa' -exec cat {} + | gzip > "${meta.alias}.transcriptome.fa.gz"
    find gffs/ -name '*.gff' -exec cat {} + |grep -v '^#' | gzip > "${meta.alias}.transcriptome.gff.gz"
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
    """
    echo packing images
    """
}


workflow process_bams {
    take:
        bam
        extracted_barcodes
        gtf
        bc_longlist_dir
        ref_genome_fasta
        ref_genome_idx
    main:
        chr_gtf = split_gtf_by_chroms(gtf)
            .flatten()
            .map {file -> 
                // create [chr, gtf]
                tuple(file.baseName, file)}

        get_contigs(bam)

        contigs = get_contigs.out.contigs
            .splitCsv().map{it -> [it[0][0], it[1]]}

        // Keep only the contigs that are referenced in the gtf
        contigs = chr_gtf
            .cross(contigs) // -> [[ chr, chr.gtf], [chr, meta]]
            // [meta, chr, chr.gtf]
            .map {chr_gtf, chr_meta -> [chr_meta[1], chr_meta[0], chr_gtf[1]]}

        generate_whitelist(
            extracted_barcodes.groupTuple(),
            bc_longlist_dir)

        assign_barcodes(
             generate_whitelist.out.whitelist
            .cross(extracted_barcodes)
            .map {it ->
                    meta = it[0][0]
                    whitelist = it[0][1]
                    barcodes = it[1][1]
                    [meta, whitelist, barcodes]
                })

        // Combine the BAM and tags chunks
        combine_bams_and_tags(
            bam.groupTuple()
                .join(assign_barcodes.out.tags.groupTuple()))

        // Split the tags by chromosome
        chr_tags = combine_bams_and_tags.out.merged_tags
            .transpose()
            .map {meta, file -> [meta, file.baseName, file]}

        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            combine_bams_and_tags.out.merged_bam
                .combine(chr_gtf))

        align_to_transcriptome(stringtie.out.read_tr_map)

        assign_features(
            align_to_transcriptome.out.read_tr_map
            .join(chr_tags, by: [0, 1]))

        cluster_umis(
            assign_features.out.feature_assigns
            // Join on [sample meta,chr]
            .join(chr_tags, by: [0, 1]))

        tag_bams(combine_bams_and_tags.out.merged_bam
             // cross by sample_id on the output of cluster_umis to return
             // [sample_id, chr, kit_name, bam, bai, tags.tsv]
            .cross(cluster_umis.out.read_tags)
            .map {it -> it.flatten()[0, 1, 2, 4, 5 ]})

        read_tags = combine_tag_files(
            cluster_umis.out.read_tags
             .map {it -> it[0, 2]}.groupTuple())

        final_read_tags = combine_final_tag_files(
            cluster_umis.out.final_read_tags.groupTuple())

        umi_gene_saturation(read_tags)

        construct_gene_expression(read_tags, 'gene')
        construct_transcript_expression(read_tags, 'transcript')

        if (params.plot_umaps == true) {
            umap_reduce_expression_matrix(
                Channel.from(1..params.umap_n_repeats)
                .combine(
                    construct_gene_expression.out.matrix_processed_tsv
                    .concat(
                        construct_transcript_expression.out.matrix_processed_tsv)))
             umaps = umap_reduce_expression_matrix.out.matrix_umap_tsv.groupTuple()
        }else{
            umaps = construct_gene_expression.out.matrix_processed_tsv
                // Make optinal file for each sample - [sample_id, OPTIONAL_FILE]
                .map {[it[0], file("$projectDir/data/OPTIONAL_FILE")]}
        }

        if (params.merge_bam) {
            combine_chrom_bams(tag_bams.out.tagged_bam
                .groupTuple())
            // [sample_id, bam]
            tagged_bams = combine_chrom_bams.out.bam_fully_tagged
        }else{
            tagged_bams = tag_bams.out.tagged_bam
            // [sample_id, bam, bai]
            .map {it -> it[0, 1, 2]}
            .groupTuple()
        }

    pack_images(
        generate_whitelist.out.kneeplot
       .concat(umi_gene_saturation.out.saturation_curve)
       .groupTuple())

    merge_transcriptome(
        assign_features.out.annotation.groupTuple()
            .join(stringtie.out.read_tr_map.groupTuple())
            .map{
                meta, ann_tr_gff, chr, tr_fa, ref_gtf, str_gff, fastq  ->
                [meta, tr_fa, ann_tr_gff]})

    // Tidy up channels prior to output. Keep only [meta, output]
    expresion_out = construct_gene_expression.out.matrix_processed_tsv
        .concat(
            construct_transcript_expression.out.matrix_processed_tsv,
            construct_gene_expression.out.matrix_counts_tsv,
            construct_transcript_expression.out.matrix_counts_tsv)
        .map {it -> it[0, 2]}.groupTuple()

    emit:
        results = umaps
            .join(umi_gene_saturation.out.saturation_curve)
            .join(final_read_tags)
            .join(expresion_out)
            .join(construct_gene_expression.out.mito_expression_tsv)
            .join(generate_whitelist.out.whitelist)
            .join(generate_whitelist.out.uncorrected_bc_counts)
            .join(generate_whitelist.out.kneeplot)
            .join(tagged_bams)
            .join(pack_images.out)
            .join(merge_transcriptome.out)
            .map{it -> it.flatten()}

        // Emit sperately for use in the report
        final_read_tags = final_read_tags
        plots = pack_images.out.collect{it -> it[1]}.collect()
        white_list = generate_whitelist.out.whitelist
        gene_expression = construct_gene_expression.out.matrix_processed_tsv.map {it -> it[0, 2]}
        transcript_expression = construct_transcript_expression.out.matrix_processed_tsv.map {it -> it[0, 2]}
        mitochondrial_expression = construct_gene_expression.out.mito_expression_tsv
        umap_matrices = umaps
}
