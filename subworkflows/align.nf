process call_paftools {
    label "singlecell"
    conda "${projectDir}/environment.yaml"
    input:
        path ref_genes_gtf
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    """
    paftools.js gff2bed -j $ref_genes_gtf > ref_genes.bed
    """
}

process get_chrom_sizes{
    input:
        path ref_genome_idx
    output:
        path 'chr_sizes', emit: ref_chrom_sizes
    """
    cut -f1,2 $ref_genome_idx | sort -V > chr_sizes
    """
}

process align_to_ref {
    label "singlecell"
    conda "${projectDir}/environment.yaml"
    input:
        tuple val(sample_id),
              path(stranded_fq)
        path ref_genes_bed
        path ref_genome_fasta
        path ref_chrom_sizes
    output:
        tuple val(sample_id), path("*sorted.bam"), emit: bam_sort
        tuple val(sample_id), path("*sorted.bam.bai"), emit: bam_sort_bai
    """
     minimap2 -ax splice -uf --MD -t $params.resources_mm2_max_threads \
      --junc-bed ${ref_genes_bed} \
      --secondary=no \
      ${ref_genome_fasta} ${stranded_fq} > tmp.sam && \
    samtools view --no-PG tmp.sam \
      -t ref_chrom_sizes -o unsort.bam;
    samtools sort --no-PG unsort.bam -o ${sample_id}_sorted.bam;
    samtools index ${sample_id}_sorted.bam
    """
}



// workflow module
workflow align {
    take:
        stranded_fq
        ref_genome_fasta
        ref_genes_gtf
        ref_genome_idx
    main:
        call_paftools(ref_genes_gtf)
        get_chrom_sizes(ref_genome_idx)
        align_to_ref(
            stranded_fq,
            call_paftools.out.ref_genes_bed,
            ref_genome_fasta,
            get_chrom_sizes.out.ref_chrom_sizes)
   
    emit:
        bam_sort = align_to_ref.out.bam_sort
        bam_sort_bai = align_to_ref.out.bam_sort_bai

}
