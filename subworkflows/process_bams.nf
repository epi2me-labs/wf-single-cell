import java.util.ArrayList;

process get_contigs {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              path('sample.bam'),
              path('sample.bam.bai')

    output:
        path("${sample_id}_contigs"),
        emit: contigs
    """
    samtools idxstats sample.bam \
        | gawk -v var=${sample_id} '/^[^*]/{print var,\$1}' \
        | gawk NF > "${sample_id}_contigs"
    """
}

process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    cpus 2
    input:
        tuple path("sort.bam"),
              path("sort.bam.bai"),
              val(meta),
              val(chrom)
        path "bc_longlist_dir"

    output:
        // TODO: Do not write bams. Write mapping of read_id to barcode
        tuple val(meta.sample_id), 
              path("*.bc_extract.sorted.tsv"),
              val(chrom),
              emit: bc_uncorr_tsv
        tuple val(meta.sample_id),
              path("*.uncorrected_bc_counts.tsv"), emit: barcode_counts

    """
    workflow-glue extract_barcode \
    sort.bam bc_longlist_dir/${meta['bc_long_list']}\
    -t $task.cpus \
    --kit ${meta['kit_name']} \
    --adapter1_suff_length $params.barcode_adapter1_suff_length \
    --min_barcode_qv $params.barcode_min_quality \
    --barcode_length ${meta['barcode_length']} \
    --umi_length ${meta['umi_length']} \
    --output_read_tags "${meta.sample_id}.bc_extract.sorted.tsv" \
    --output_barcode_counts "${meta.sample_id}.${chrom}.uncorrected_bc_counts.tsv" \
    --contig ${chrom}
    """
}

process generate_whitelist{
    label "singlecell"
    cpus 1
    input:
        tuple path("counts"),
              val(meta)
    output:
        tuple val(meta.sample_id), 
              path("*whitelist.tsv"), 
              emit: whitelist
        tuple val(meta.sample_id), 
              path("*kneeplot.png"), 
              emit: kneeplot
    """
    workflow-glue knee_plot \
        counts \
        --exp_cells ${meta['exp_cells']} \
        --output_whitelist "${meta.sample_id}.whitelist.tsv" \
        --output_plot "${meta.sample_id}.kneeplot.png"
    """
}

process assign_barcodes{
    label "singlecell"
    cpus 1
    input:
         tuple val(sample_id),
               path("whitelist.tsv"),
               path("extract_barcodes.tsv"),
               val(chr)
    output:
        tuple val(sample_id),
              val(chr),
              path("bc_assign_counts.tsv"),
              emit: chrom_assigned_barcode_counts
        tuple val(sample_id),
              val(chr),
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
    input:
        tuple val(sample_id),
              val(chr),
              path("chrom_feature_assigns.tsv"),
              path("read_tags.tsv")
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}_${chr}.read_tags.tsv"),
              emit: read_tags  // For BAM tagging
        tuple val(sample_id),
              path("${sample_id}_${chr}.final_tags.tsv"),
              emit: final_read_tags  // For user output
    """
    workflow-glue cluster_umis \
        --chrom ${chr} \
        --cell_gene_max_reads $params.umi_cell_gene_max_reads \
        --feature_assigns chrom_feature_assigns.tsv \
        --read_tags read_tags.tsv \
        --output_read_tags "${sample_id}_${chr}.read_tags.tsv" \
        --workflow_output "${sample_id}_${chr}.final_tags.tsv"
    """
}

process tag_bams {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path("align.bam"),
              path("align.bam.bai"),
              path('tags.tsv')
    output:
         tuple val(sample_id),
              path("${sample_id}.${chr}.tagged.bam"),
              path("${sample_id}.${chr}.tagged.bam.bai"),
              emit: tagged_bam
    """
    workflow-glue tag_bam \
        --in_bam align.bam \
        --tags tags.tsv \
        --out_bam ${sample_id}.${chr}.tagged.bam

    samtools index ${sample_id}.${chr}.tagged.bam
    """
}

process combine_tag_files {
    // A file named 'read_tags' is a required output.
    // collectFile's 'name' argument does not work when being applied to
    // a channel that returns tuples. It groups and names according to the
    // first element of the tuple. Hence this process.
    label "singlecell"
    input:
        tuple val(sample_id),
              path("tags*.tsv")
    output:
        tuple val(sample_id),
              path("${sample_id}_read_tags.tsv")
    """
    awk 'FNR>1 || NR==1' *.tsv > "${sample_id}_read_tags.tsv"
    """
}

process combine_final_tag_files {
    // Combine the final
    label "singlecell"
    input:
        tuple val(sample_id),
              path("tags*.tsv")
    output:
        tuple val(sample_id),
              path("${sample_id}.read_tags.tsv")
    """
    awk 'FNR>1 || NR==1' *.tsv > "${sample_id}.read_tags.tsv"
    """
}

process combine_uncorrect_bcs {
    label "singlecell"
    input:
        tuple val(sample_id),
              path("uncorrected_bcs*.tsv")
        output:
            tuple val(sample_id),
                  path("${sample_id}.uncorrected_bc_counts.tsv")
    """
    cat *.tsv > "${sample_id}.uncorrected_bc_counts.tsv"
    """
}



process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "singlecell"
    cpus Math.min(8, params.max_threads)
    input:
        tuple val(sample_id), 
              path(chrom_bams),
              path('chrom.bam.bai')
    output:
        tuple val(sample_id), 
              path("*tagged.sorted.bam"), 
              path("*tagged.sorted.bam.bai"),
              emit: bam_fully_tagged
    """
    samtools merge -@ ${task.cpus} -o "${sample_id}.tagged.sorted.bam" ${chrom_bams}; 
    samtools index -@ ${task.cpus} "${sample_id}.tagged.sorted.bam";
    """
}


process stringtie {
    label "singlecell"
    cpus Math.max(params.max_threads / 4, 4.0)
    input:
        path 'ref_genome.fa'
        path 'ref_genome.fa.fai'
        tuple val(sample_id),
              path("align.bam"),
              path("align.bam.bai"),
              val(chr),
              path("chr.gtf")
    output:
        tuple val(sample_id),
              val(chr),
              path("transcriptome.fa"),
              path("chr.gtf"),
              path("stringtie.gff"),
              path("reads.fastq"),
              emit: read_tr_map
    """
    # Build transcriptome. 
    workflow-glue process_bam_for_stringtie align.bam ${chr}  \
        | tee >(stringtie -L  -p ${task.cpus} -G chr.gtf -l stringtie \
            -o stringtie.gff - ) \
        | samtools fastq > reads.fastq

    # Get transcriptome sequence
    gffread -g ref_genome.fa -w transcriptome.fa stringtie.gff
    """
}


process align_to_transcriptome {
    label "singlecell"
    cpus Math.max(params.max_threads / 2, 4.0)
    input:
        tuple val(sample_id),
              val(chr),
              path('transcriptome.fa'),
              path('chr.gtf'),
              path('stringtie.gff'),
              path("reads.fastq")
    output:
        tuple val(sample_id),
              val(chr),
              path("chr.gtf"),
              path("read_query_tr_map.tsv"),
              path('stringtie.gff'),
              emit: read_tr_map
    """
    echo "read_id\tqry_id\n" > read_query_tr_map.tsv;
    minimap2 -ax map-ont \
        --end-bonus 10 \
        -t $task.cpus \
        transcriptome.fa \
        reads.fastq \
        | samtools view -F 2304 \
        |  gawk 'BEGIN{OFS="\t";} /^[^@]/ {print \$1,\$3}' \
        >> read_query_tr_map.tsv;
    """
}


process assign_features {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path("chr.gtf"),
              path("query_transcript_read_assign.tsv"),
              path('stringtie.gff'),
              path('tags.tsv')
    output:
        tuple val(sample_id),
              val(chr),
              path("${sample_id}.${chr}.feature_assigns.tsv"),
              emit: feature_assigns
    """
    if [ \$(wc -l <read_feature_map.tsv) -eq 0 ]; then
        touch ${sample_id}.${chr}_empty.transcript_assigns.tsv  
    else
        gffcompare -o gffcompare -r chr.gtf stringtie.gff
        
        workflow-glue assign_features \
            --query_transcript_read_assign query_transcript_read_assign.tsv \
            --gffcompare_tmap gffcompare.stringtie.gff.tmap \
            --gtf chr.gtf \
            --tags tags.tsv \
            --output "${sample_id}.${chr}.feature_assigns.tsv" \
            --min_mapq ${params.gene_assigns_minqv}
    fi
    """
}

process umi_gene_saturation {
    label "singlecell"
    cpus 4
    input:
        tuple val(sample_id),
              path("read_tags.tsv")
    output:
        tuple val(sample_id),
              path("*saturation_curves.png"),
              emit: saturation_curve
    """
    workflow-glue calc_saturation \
        --output "${sample_id}.saturation_curves.png" \
        --threads ${task.cpus} \
        --read_tags read_tags.tsv
    """
}

process construct_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("read_tags.tsv")
    output:
        tuple val(sample_id), 
              path("*gene_expression.counts.tsv"),
              path("*transcript_expression.counts.tsv"),
              emit: matrix_counts_tsv
    """
    workflow-glue gene_expression \
        --output_prefix "${sample_id}" \
        --read_tags read_tags.tsv
    """
}

process process_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("gene_matrix_counts.tsv"),
              path("transcript_matrix_counts.tsv")
    output:
        tuple val(sample_id),
              val('gene'),
              path("${sample_id}.gene_expression.processed.tsv"),
              emit: gene_matrix_processed_tsv
        tuple val(sample_id),
              val('transcript'),
              path("${sample_id}.transcript_expression.processed.tsv"),
              emit: transcript_matrix_processed_tsv
        tuple val(sample_id),
              path("${sample_id}.gene_expression.mito.tsv"),
              emit: mito_expression_tsv
    """
    workflow-glue process_matrix \
        --min_genes $params.matrix_min_genes \
        --min_cells $params.matrix_min_cells \
        --max_mito $params.matrix_max_mito \
        --mito_prefix ${params.mito_prefix} \
        --norm_count $params.matrix_norm_count \
        --gene_counts gene_matrix_counts.tsv \
        --sample_id ${sample_id} \
        --transcript_counts transcript_matrix_counts.tsv
    """
}


process umap_reduce_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(repeat_num),
              val(sample_id),
              val(data_type),
              path(matrix)
    output:
        tuple val(sample_id),
                path("${data_type}_umap_${repeat_num}.tsv"),
                emit: matrix_umap_tsv
    """
    workflow-glue umap_reduce \
        --output ${data_type}_umap_${repeat_num}.tsv \
        ${matrix}
    """
}


process pack_images {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("images_${sample_id}/*")
    output:
         tuple val(sample_id),
              path("images_${sample_id}")
    """
    echo packing images
    """
}


workflow process_bams {
    take:
        bam
        meta
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
            .splitCsv(sep: " ")
        
        // Keep only the contigs that are referenced in the gtf
        contigs = chr_gtf.map {[it[1], it[0]]}
            .cross(contigs) {it -> it[1]}
            .map {it -> it.flatten()[2, 3]}

        extract_barcodes(
            bam
            .cross(
                meta
                .cross(contigs).map{it -> it.flatten()})
                .map{it -> it.flatten()[1, 2, 4, 6]},
            bc_longlist_dir)

        un_corr_bcs = combine_uncorrect_bcs(extract_barcodes.out.barcode_counts)

        generate_whitelist(
            extract_barcodes.out.barcode_counts
            .collectFile()
            .map {it -> tuple(it.getSimpleName(), it)}
            .join(meta).map {it -> it.tail()}) // Remove sample_id


       assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extract_barcodes.out.bc_uncorr_tsv)
            .map {it -> it.flatten()[0, 1, 3, 4]})

        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            bam.combine(chr_gtf))
        
        align_to_transcriptome(stringtie.out.read_tr_map)

        assign_features(
            align_to_transcriptome.out.read_tr_map
            .join(assign_barcodes.out.tags, by: [0, 1]))

        cluster_umis(
            assign_features.out.feature_assigns
            // Join on sample_id,chr
            .join(assign_barcodes.out.tags, by: [0, 1]))

        tag_bams(
            bam.cross(
                cluster_umis.out.read_tags
            ).map {it ->it.flatten()[0, 4, 1, 2, 5]})

        read_tags = combine_tag_files(
            cluster_umis.out.read_tags
             .map {it -> it[0, 2]}.groupTuple())

       final_read_tags = combine_final_tag_files(
            cluster_umis.out.final_read_tags.groupTuple())

        umi_gene_saturation(read_tags)

        construct_expression_matrix(read_tags)

        process_expression_matrix(
            construct_expression_matrix.out.matrix_counts_tsv)

        if (params.plot_umaps == true) {
            umap_reduce_expression_matrix(
                Channel.from(1..params.umap_n_repeats)
                .combine(
                    process_expression_matrix.out.gene_matrix_processed_tsv
                    .concat(
                        process_expression_matrix.out.transcript_matrix_processed_tsv)))
             umaps = umap_reduce_expression_matrix.out.matrix_umap_tsv.groupTuple()
        }else{
            umaps = process_expression_matrix.out.gene_matrix_processed_tsv
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

    // Tidy up channels prior to output
    proc_expresion_out = process_expression_matrix.out.gene_matrix_processed_tsv
        .concat(process_expression_matrix.out.transcript_matrix_processed_tsv)
        .map {it -> it[0, 2]}.groupTuple()

    emit:
        results = umaps
            .mix(umi_gene_saturation.out.saturation_curve)
            .join(final_read_tags)
            .join(construct_expression_matrix.out)
            .join(proc_expresion_out)
            .join(process_expression_matrix.out.mito_expression_tsv)
            .join(generate_whitelist.out.whitelist)
            .join(generate_whitelist.out.kneeplot)
            .join(tagged_bams)
            .join(combine_uncorrect_bcs.out)
            .join(pack_images.out)
            .map{it -> it.flatten()}
        
        // Emit sperately for use in the report
        final_read_tags = final_read_tags
        plots = pack_images.out.collect{it -> it[1]}.collect()
        white_list = generate_whitelist.out.whitelist
        gene_expression = process_expression_matrix.out.gene_matrix_processed_tsv.map {it -> it[0, 2]}
        transcript_expression = process_expression_matrix.out.transcript_matrix_processed_tsv.map {it -> it[0, 2]}
        mitochondrial_expression = process_expression_matrix.out.mito_expression_tsv
        umap_matrices = umaps
}
