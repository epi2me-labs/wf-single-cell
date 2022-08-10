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
    "cut -f1,2 $ref_genome_idx | sort -V > chr_sizes"
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
        path "${sample_id}_sorted.bam", emit: BAM_SORT
        path "${sample_id}_sorted.bam.bai", emit: BAM_SORT_BAI
    """
    minimap2 -ax splice -uf --MD -t 2 \
      --junc-bed ${ref_genes_bed} \ 
      --secondary=no \ 
      ${REF_GENOME_FASTA} ${STRANDED_FQ} > \
      tmp.sam && \
    samtools view --no-PG tmp.sam \
      -t REF_CHROM_SIZES -o unsort.bam;
    samtools sort --no-PG ${sample_id}_unsort.bam -o ${sample_id}_sorted.bam;
    samtools index ${sample_id}_sorted.bam.bai
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
        d = {it ->
        /* Harmonize tuples
        output:
            tuple val(sample_id), path('*.gff')
        When there are multiple paths, will emit:
            [sample_id, [path, path ..]]
        when there's a single path, this:
            [sample_id, path]
        This closure makes both cases:
            [[sample_id, path][sample_id, path]].
        */
            if (it[1].getClass() != java.util.ArrayList){
                // If only one path, `it` will be [sample_id, path]
                return [it]
            }
            l = [];
            for (x in it[1]){
                l.add(tuple(it[0], x))
            }
            return l
        }
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
