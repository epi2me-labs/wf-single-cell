process call_paftools {
    label "singlecell"
    memory "2 GB"
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
    memory "1 GB"
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
    cpus params.threads
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
    cpus params.threads
    memory "32 GB"
    input:
        tuple val(group_index),
              val(meta),
              path("reads.fastq")
        path "genome_index.mmi"
        path "ref_genes.bed"
        path "ref_chrom_sizes.tsv"
    output:
        tuple val(meta),
              path("sorted.bam"), 
              path("sorted.bam.bai"), 
              emit: bam_sort
        tuple val(group_index),
              val(meta),
              path("bam_info.tsv"),
              emit: bam_info
    script:
    def view_threads = 1
    def sort_threads = 3
    def mm2_threads = Math.max(task.cpus - view_threads - sort_threads, 4)
    """
    minimap2 -ax splice -uf --secondary=no --MD -t $mm2_threads \
        --junc-bed ref_genes.bed  \
        --cap-kalloc 100m --cap-sw-mem 50m \
        genome_index.mmi reads.fastq \
    | samtools view -@ $view_threads -b --no-PG -t ref_chrom_sizes - \
    | tee >(samtools sort -@ $sort_threads --write-index -o "sorted.bam"##idx##"sorted.bam.bai" --no-PG  -) \
    | seqkit bam -F - 2> bam_info.tsv
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
