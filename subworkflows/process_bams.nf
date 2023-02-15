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

process assign_genes {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path("chrom_bc.bed"),
              path("chrom.gtf")
    output:
        tuple val(sample_id),
              val(chr),
              path("*.read.gene_assigns.tsv"),
              emit: chrom_tsv_gene_assigns
    """
    workflow-glue assign_genes \
    --output "${sample_id}_${chr}.read.gene_assigns.tsv" \
    chrom_bc.bed chrom.gtf
    """
}

process cluster_umis {
    label "singlecell"
    cpus params.umi_cluster_max_threads
    input:
        tuple val(sample_id),
              path("align.bam"),
              path("align.bam.bai"),
              val(chr),
              path("chrom_gene_assigns.tsv"),
              path("chrom_tr_assigns.tsv"),
              path(bc_ur_tags)
    output:
         tuple val(sample_id),
              path("*.tagged.bam"),
              path("*.tagged.bam.bai"),
              emit: bam_bc_bai
        tuple val(sample_id),
              path("*tags.tsv"),
              emit: tags
    """
    workflow-glue cluster_umis \
    align.bam \
    --threads $task.cpus \
    --chrom ${chr} \
    --ref_interval $params.umi_genomic_interval \
    --cell_gene_max_reads $params.umi_cell_gene_max_reads \
    --gene_assigns chrom_gene_assigns.tsv \
    --transcript_assigns chrom_tr_assigns.tsv \
    --bc_ur_tags ${bc_ur_tags} \
    --output "${sample_id}_${chr}.tagged.bam" \
    --output_read_tags "${sample_id}_${chr}.read_tags.tsv" \
    --sample_id ${sample_id}

    samtools index "${sample_id}_${chr}.tagged.bam"
    """
}

process combine_tag_files {
    // collectFile does 'name' argument does not work when being applied to
    // A channel that returns tuples. It groups and names according to the
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
              path("read_order.tsv"),
              path("reads.fastq"),
              emit: read_tr_map
    """
    # Build transcriptome. 
    workflow-glue process_bam_for_stringtie align.bam ${chr}  \
        | tee >(stringtie -L  -p ${task.cpus} -G chr.gtf -l stringtie \
            -o stringtie.gff - ) \
        | samtools fastq - \
        | tee reads.fastq \
        | seqkit seq -n > read_order.tsv 
    
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
              path('read_order.tsv'),
              path("reads.fastq")
    output:
        tuple val(sample_id),
              val(chr),
              path("chr.gtf"),
              path("*read_query_tr_map.tsv"), 
              path('stringtie.gff'),
              path("read_order.tsv"),
              emit: read_tr_map
    """  
    minimap2 -ax map-ont \
        --end-bonus 10 \
        -t $task.cpus \
        transcriptome.fa \
        reads.fastq \
        | samtools view -F 2304 \
        |  gawk 'BEGIN{OFS="\t";} /^[^@]/ {print \$1,\$3}' \
        > read_query_tr_map.tsv;
    """
}


process assign_transcripts {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path("chr.gtf"),
              path("read_query_tr_map.tsv"),
              path('stringtie.gff'),
              path('read_order.tsv')
    output:
        tuple val(sample_id),
              val(chr),
              path("*transcript_assigns.tsv"), 
              emit: transcript_assigns
    """
    if [ \$(wc -l <read_query_tr_map.tsv) -eq 0 ]; then
        touch ${sample_id}.${chr}_empty.transcript_assigns.tsv  
    else
        gffcompare -o gffcompare -r chr.gtf stringtie.gff
        
        workflow-glue isoform_read_mapping \
            --read_tr_map read_query_tr_map.tsv \
            --all_read_ids read_order.tsv \
            --gffcompare_tmap gffcompare.stringtie.gff.tmap \
            --output "${sample_id}.${chr}.transcript_assigns.tsv"
    fi
    """
}

process umi_gene_saturation {
    label "singlecell"
    cpus 1
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
        read_tags.tsv
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
              path("*gene_expression.processed.tsv"),
              path("*transcript_expression.processed.tsv"),
              emit: matrix_processed_tsv
        tuple val(sample_id),
              path("*gene_expression.mito.tsv"),
              emit: matrix_mito_tsv
    """
    workflow-glue process_matrix \
    --min_genes $params.matrix_min_genes \
    --min_cells $params.matrix_min_cells \
    --max_mito $params.matrix_max_mito \
    --mito_prefix ${params.mito_prefix} \
    --norm_count $params.matrix_norm_count \
    --output_prefix ${sample_id} \
    --gene_counts gene_matrix_counts.tsv \
    --transcript_counts transcript_matrix_counts.tsv
    """
}


process umap_reduce_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("gene_matrix_processed.tsv"),
              path("transcript_matrix_processed.tsv")
    output:
         tuple val(sample_id),
              path("*gene_expression*umap.tsv"),
              path("*transcript_expression*umap.tsv"),
            //   path("gene_matrix_processed.tsv"),
            //   path("transcript_matrix_processed.tsv"),
              emit: matrix_umap_tsv
    """
    workflow-glue umap_reduce \
        --output_prefix "${sample_id}.gene_expression" \
        gene_matrix_processed.tsv

    workflow-glue umap_reduce \
        --output_prefix "${sample_id}.transcript_expression" \
        transcript_matrix_processed.tsv
    """
}


process umap_plot_total_umis {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(gene_matrix_umap),
              path(transcript_matrix_umap),
              path("gene_matrix_processed.tsv"),
              path("transcript_matrix_processed.tsv")
    output:
          tuple val(sample_id),
              path("*.genes*.png"), 
              emit: gene_umap_plot_total
          tuple val(sample_id),
              path("*transcripts*.png"), 
              emit: transcript_umap_plot_total

    """
    workflow-glue plot_umap \
        --output_prefix "${sample_id}.umap.genes.total" \
        --umap ${gene_matrix_umap} \
        --full_matrix gene_matrix_processed.tsv
    
    workflow-glue plot_umap \
        --output_prefix "${sample_id}.umap.transcripts.total" \
        --umap ${transcript_matrix_umap} \
        --full_matrix transcript_matrix_processed.tsv
    """
}


process annotate_umap_genes {
    // TODO: make a channle of input genes for thes process
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_umap_gene),
              path(matrix_umap_transcript),
              path("matrix_processed_gene.tsv"),
              path("matrix_processed_transcript.tsv"),
              val(gene)
    output:
        tuple val(sample_id),
              path("*.png"),
              emit: umaps
    script:
    """
    workflow-glue plot_umap \
        --gene $gene \
        --output_prefix "${sample_id}.umap.gene_annotate.${gene}" \
        --umap ${matrix_umap_gene} \
        --full_matrix matrix_processed_gene.tsv
    """
}

process pack_images {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path("images_${sample_id}/*")
    output:
        path "images_${sample_id}"
    """
    echo packing images
    """
}

process umap_plot_mito_genes {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(gene_matrix_umaps),
              path("matrix_mito_processed.tsv")
    output:
        tuple val(sample_id),
              path("*.png"), 
              emit: matrix_umap_plot_mito
    """
    workflow-glue plot_umap \
        --mito_genes \
        --output_prefix "${sample_id}.umap.mitochondrial" \
        --umap ${gene_matrix_umaps} \
        --full_matrix matrix_mito_processed.tsv
    """
}
    
workflow process_bams {
    take:
        bam
        chr_beds
        meta
        gtf
        umap_genes
        bc_longlist_dir
        ref_genome_fasta
        ref_genome_idx
   
    main:
        beds = chr_beds
            .flatMap{
                l = []
                for (i=0; i<it[1].size(); i++){
                    sample_id = it[0]
                    chr = it[1][i].getSimpleName()
                    bed_file = it[1][i]
                    l.add(tuple(chr, sample_id, bed_file))
                }
                return l
            }
        
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
        
        // barcodes()
        extract_barcodes(
            bam
            .cross(
                meta
                .cross(contigs).map{it -> it.flatten()})
                .map{it -> it.flatten()[1, 2, 4, 6]},
            bc_longlist_dir)

        generate_whitelist(
            extract_barcodes.out.barcode_counts
            .collectFile()
            .map {it -> tuple(it.getSimpleName(), it)}
            .join(meta).map {it -> it.tail()} // Remove sample_id
        )

        assign_barcodes(
            generate_whitelist.out.whitelist
            .cross(extract_barcodes.out.bc_uncorr_tsv)
            .map {it -> it.flatten()[0, 1, 3, 4]}
       )

        // combine all chr bams with chr gtfs
        chr_beds_gtf = chr_gtf.cross(
            beds)
            // [sample_id, chr, bed, gtf]
            .map {it -> it.flatten()[3, 0, 4, 1]}
            

        assign_genes(chr_beds_gtf)

        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            bam.combine(chr_gtf))
        
        align_to_transcriptome(
            stringtie.out.read_tr_map
        )
        
        assign_transcripts {
            align_to_transcriptome.out.read_tr_map
        }

        cluster_umis(
            bam.cross(
            // Cross on sample_id
            assign_genes.out.chrom_tsv_gene_assigns
            //join on sample_id + chr
            .join(assign_transcripts.out.transcript_assigns, by: [0, 1])
            .join(assign_barcodes.out.tags, by: [0, 1]))
            .map {it -> it.flatten()[0, 1, 2, 4, 5, 6, 7]})

        read_tags = combine_tag_files(
            cluster_umis.out.tags.groupTuple())

        umi_gene_saturation(read_tags)

        construct_expression_matrix(read_tags)

        process_expression_matrix(
            construct_expression_matrix.out.matrix_counts_tsv)

         umap_reduce_expression_matrix(
            process_expression_matrix.out.matrix_processed_tsv)

         umap_plot_total_umis(
            umap_reduce_expression_matrix.out.matrix_umap_tsv
            .join(process_expression_matrix.out.matrix_processed_tsv))

         genes_to_plot = Channel.fromPath(umap_genes)
             .splitCsv()
        
         annotate_umap_genes(
             umap_reduce_expression_matrix.out.matrix_umap_tsv
             .join(process_expression_matrix.out.matrix_processed_tsv)
             // Combine each umap with a gene to annotate with
             .transpose()
             .combine(genes_to_plot)
             .transpose())
        
        umap_plot_mito_genes(
            umap_reduce_expression_matrix.out.matrix_umap_tsv
            .join(process_expression_matrix.out.matrix_mito_tsv)
            .map {it -> it[0, 1, 3]})

        
        if (params.merge_bam) {
            combine_chrom_bams(cluster_umis.out.bam_bc_bai
                .groupTuple())
            // [sample_id, bam]
            tagged_bams = combine_chrom_bams.out.bam_fully_tagged

        }else{
            tagged_bams = cluster_umis.out.bam_bc_bai
            // [sample_id, bam, bai]
            .map {it -> it[0, 1, 2]}
            .groupTuple()
        }

    // Select the first replicate of each umap for report plotting.
    umaps = umap_plot_total_umis.out.gene_umap_plot_total
        .concat(
            umap_plot_total_umis.out.transcript_umap_plot_total,
            umap_plot_mito_genes.out,
            annotate_umap_genes.out.umaps.groupTuple())
        .flatMap {it -> 
                l = []
                for (i=0; i<it[1].size(); i++){
                    l.add([it[0], it[1][i]])
                }
                return l
        }.filter {it -> it[1] ==~ /.*_0\.png$/}


    pack_images(
           umaps
            .concat(generate_whitelist.out.kneeplot,
                umi_gene_saturation.out.saturation_curve)
            .groupTuple())

     emit:
        results = umi_gene_saturation.out.saturation_curve
             .join(read_tags)
             .join(construct_expression_matrix.out)
             .join(process_expression_matrix.out.matrix_processed_tsv)
             .join(process_expression_matrix.out.matrix_mito_tsv)
             .join(generate_whitelist.out.whitelist)
             .join(generate_whitelist.out.kneeplot)
             .join(tagged_bams)
             .join(extract_barcodes.out.barcode_counts)
             .join(annotate_umap_genes.out.umaps)
             .join(umap_reduce_expression_matrix.out.matrix_umap_tsv)
             .join(umap_plot_total_umis.out.gene_umap_plot_total)
             .join(umap_plot_mito_genes.out.matrix_umap_plot_mito)
             .join(umap_plot_total_umis.out.transcript_umap_plot_total)
             .map{it -> it.flatten()}
        
        // Emit these sperately for use in the report
        read_tags = read_tags
        plots = pack_images.out.collect()
        white_list = generate_whitelist.out.whitelist
}
