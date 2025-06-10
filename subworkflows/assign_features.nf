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


process stringtie {
    label "singlecell"
    cpus params.threads
    // Memory usage for this process is usually less than 3GB, but some cases it may go over this.
    memory { 3.GB * task.attempt }
    maxRetries 3
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
    memory { 1.0.GB.toBytes() + (tags.size() * 4 ) }
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
	touch test
    workflow-glue assign_features \
        tr_align.bam \
        gffcompare.stringtie.gff.tmap \
        chr.gtf \
        tags.tsv \
        feature_assigns.tsv \
        --min_mapq ${params.gene_assigns_minqv}
    """
}



workflow assign_features_with_stringtie {
    take:
        bam
        chr_tags
        ref_gtf
        ref_genome_fasta
        ref_genome_idx
    main:
         chr_gtf = split_gtf_by_chroms(ref_gtf)
            .flatten()
            .map {fname -> tuple(fname.baseName, fname)} // [chr, gtf]
        
        stringtie(
            ref_genome_fasta,
            ref_genome_idx,
            bam.combine(chr_gtf))

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

    emit:
        feaure_assignmnets = assign_features.out.feature_assigns
        annotation = assign_features.out.annotation
        read_to_transcript_map = stringtie.out.read_tr_map

}
