process call_paftools {
    label "singlecell"
    cpus 1
    input:
        path "ref_genes.gtf"
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    """
    paftools.js gff2bed -j ref_genes.gtf > ref_genes.bed
    """
}

process get_chrom_sizes{
    label "singlecell"
    cpus 1
    input:
        path "ref_genome.fai"
    output:
        path 'chr_sizes', emit: ref_chrom_sizes
    """
    cut -f1,2 ref_genome.fai | sort -V > chr_sizes
    """
}

process build_minimap_index{
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    cpus params.resources_mm2_max_threads
    memory '16 GB'
    input:
        path "reference.fa"
    output:
        path "genome_index.mmi", emit: index
    script:
    """
    minimap2 -t ${task.cpus} -I 16G -d "genome_index.mmi" "reference.fa"
    """
}

process align_to_ref {
    label "singlecell"
    cpus params.resources_mm2_max_threads
    input:
        tuple val(meta),
              val(chunk_id),
              path("reads.fastq")
        path "genome_index.mmi"
        path "ref_genes.bed"
        path "ref_chrom_sizes.tsv"
    output:
        tuple val(meta),
              path("${meta.alias}_sorted.bam"), 
              path("${meta.alias}_sorted.bam.bai"), 
              emit: bam_sort
        tuple val(meta),
              val(chunk_id),
              path("bam_info.tsv"),
              emit: bam_info
    """
     minimap2 -ax splice -uf --secondary=no --MD -t $task.cpus \
      --junc-bed ref_genes.bed $params.resources_mm2_flags  \
      genome_index.mmi reads.fastq \
        | samtools view -b --no-PG -t ref_chrom_sizes - \
        | tee >(samtools sort -@ 2 --no-PG  - > "${meta.alias}_sorted.bam") \
        | seqkit bam -F - 2> bam_info.tsv

    samtools index -@ ${task.cpus} "${meta.alias}_sorted.bam"
    """
}


// workflow module
workflow align {
    take:
        stranded_fq
        ref_genome
        ref_genome_idx
        ref_genes_gtf
    main:
        call_paftools(ref_genes_gtf)
        get_chrom_sizes(ref_genome_idx)
        build_minimap_index(ref_genome)
        align_to_ref(
            stranded_fq,
            build_minimap_index.out.index,
            call_paftools.out.ref_genes_bed,
            get_chrom_sizes.out.ref_chrom_sizes)
    emit:
        bam_sort = align_to_ref.out.bam_sort
        bam_info = align_to_ref.out.bam_info
}
