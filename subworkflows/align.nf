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

process align_to_ref {
    label "singlecell"
    cpus params.resources_mm2_max_threads
    input:
        tuple val(sample_id),
              path("reads.fastq")
        path "ref_genome.fasta"
        path "ref_genes.bed"
        path "ref_chrom_sizes.tsv"
    output:
        tuple val(sample_id), 
            path("*sorted.bam"), 
            path("*sorted.bam.bai"), 
            emit: bam_sort
        tuple val(sample_id),
              path("beds/*.bed"),
              emit: chr_beds
    """
     minimap2 -ax splice -uf --MD -t $task.cpus \
      --junc-bed ref_genes.bed $params.resources_mm2_flags  \
      ref_genome.fasta reads.fastq* \
        | samtools view -F 2304 -b --no-PG -t ref_chrom_sizes - \
        | samtools sort -@ 2 --no-PG  - \
            | tee "${sample_id}_sorted.bam" \
        | bedtools bamtobed -i stdin \
        | gawk '/^[^#]/ {print>\$1".bed"}'
    samtools index -@ ${task.cpus} "${sample_id}_sorted.bam"

    mkdir beds
    mv *.bed beds
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
        align_to_ref(
            stranded_fq.groupTuple(),
            ref_genome,
            call_paftools.out.ref_genes_bed,
            get_chrom_sizes.out.ref_chrom_sizes)
    emit:
        bam_sort = align_to_ref.out.bam_sort
        chr_beds =  align_to_ref.out.chr_beds
}
