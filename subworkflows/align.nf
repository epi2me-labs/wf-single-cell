process call_paftools {
    input:
        path REF_GENES_GTF
    output:
        path "ref_genes.bed", emit: REF_GENES_BED
    """
    paftools.js gff2bed -j $REF_GENES_GTF > ref_genes.bed
    """
}

process get_chrom_sizes{
    input:
        // genome=str(REF_GENOME_FASTA) + ".fai",
        path ref_genome_idx
    output:
        // chrsizes=REF_CHROM_SIZES,
        path 'chr_sizes', emit: REF_CHROM_SIZES
    """
    cut -f1,2 $ref_genome_idx | sort -V > chr_sizes
    """
}


// def get_split_ont_align_mem_gb(wildcards, threads):
//     return config["RESOURCES_MM2_MEM_GB"] / threads


process align_to_ref {
    input:
        tuple val(sample_id),
              path(STRANDED_FQ)
        path ref_genes_bed
        path REF_GENOME_FASTA
        path REF_CHROM_SIZES
    output:
        tuple val(sample_id), path("*sorted.bam"), emit: BAM_SORT
        tuple val(sample_id), path("*sorted.bam.bai"), emit: BAM_SORT_BAI
    """
     minimap2 -ax splice -uf --MD -t $params.RESOURCES_MM2_MAX_THREADS \
      --junc-bed ${ref_genes_bed} \
      --secondary=no \
      ${REF_GENOME_FASTA} ${STRANDED_FQ} > tmp.sam && \
    samtools view --no-PG tmp.sam \
      -t REF_CHROM_SIZES -o unsort.bam;
    samtools sort --no-PG unsort.bam -o ${sample_id}_sorted.bam;
    samtools index ${sample_id}_sorted.bam
    """
}



// workflow module
workflow align {
    take:
        STRANDED_FQ
        REF_GENOME_FASTA
        REF_GENES_GTF
        ref_genome_idx
    main:
        call_paftools(REF_GENES_GTF)
        get_chrom_sizes(ref_genome_idx)
        align_to_ref(
            STRANDED_FQ,
            call_paftools.out.REF_GENES_BED,
            REF_GENOME_FASTA,
            get_chrom_sizes.out.REF_CHROM_SIZES)
   
    emit:
        BAM_SORT = align_to_ref.out.BAM_SORT
        BAM_SORT_BAI = align_to_ref.out.BAM_SORT_BAI

}
