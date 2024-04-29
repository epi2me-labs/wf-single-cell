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


process get_chrom_sizes {
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


process build_minimap_index {
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


process call_adapter_scan {
    label "singlecell"
    cpus params.threads
    // memory here is taken by minimap2. Having merged the three steps into one,
    // we have have prehps reduced parallelism in the workflow because in some setups
    // it might be the case that multiple tasks of the first two steps cannot now run
    // in parallel. The resolution to that would be to make the first two steps do
    // better parallelism. The advantage here is not having to write to disk, stage files
    // and read from disk between the steps (creating a lot of big temporary files).
    // 
    // peak RSS for aligning this data is robustly <12.4 GB with human reference. Set
    // a little more and do a retry
    memory {15.GB * task.attempt}
    maxRetries 1
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta), path(chunk, stageAs: 'chunk.fq.gz')
        path "bc_longlist_dir"
        path "genome_index.mmi"
        path "ref_genes.bed"
        path "ref_chrom_sizes.tsv"
    output:
        tuple val(meta), path("adapt_scan.tsv"), emit: read_config_chunked
        tuple val(meta), path("read_tags.tsv"), emit: read_tags
        tuple val(meta), path("high_quality_bc_counts.tsv"), emit: barcode_counts
        tuple val(meta), path("sorted.bam"), path("sorted.bam.bai"), emit: bam_sort
    script:
    def fl = params.full_length_only ? "--keep_fl_only": ""
    // alignment is the real bottleneck here, don't worry about threads
    // for sorting. Just subtract 1 thread as a loose bookeeping. Note the
    // hidden call to vsearch in the first program: the pipe doesn't get
    // going until thats finished. vsearch appears to use all the juice
    // it can squeeze.
    // We set -K option to minimap2 as default appears too large to
    // stream data effectively (it just blocks with defaults waiting for
    // more data). The effectiveness of this is not clear.
    def mm2_threads = task.cpus - 1
    """
    export POLARS_MAX_THREADS=$task.cpus

    workflow-glue adapter_scan_vsearch \
        chunk.fq.gz \
        --kit ${meta['kit_name']} \
        --output_tsv  "adapt_scan.tsv" \
        ${fl} \
    | \
    workflow-glue extract_barcode \
        - \
        bc_longlist_dir/${meta['bc_long_list']} \
        --kit ${meta["kit_name"]} \
        --adapter1_suff_length $params.barcode_adapter1_suff_length \
        --min_barcode_qv $params.barcode_min_quality \
        --barcode_length ${meta['barcode_length']} \
        --umi_length ${meta['umi_length']} \
        --output_read_tags "bc_extract.tsv" \
        --output_barcode_counts "high_quality_bc_counts.tsv" \
    | \
    minimap2 -ax splice -uf --secondary=no --MD \
        -t $mm2_threads -K 10M \
        --junc-bed ref_genes.bed  \
        --cap-kalloc 100m --cap-sw-mem 50m \
        genome_index.mmi - \
    | samtools view -u --no-PG -t ref_chrom_sizes - \
    | tee >(samtools sort --write-index -o "sorted.bam"##idx##"sorted.bam.bai" --no-PG  -) \
    | seqkit bam -F - 2> bam_info.tsv

    # TODO: improve this with pipes?
    csvtk cut -tlf Read,Pos,EndPos,Ref,MapQual bam_info.tsv > bam_info_cut.tsv
    # Left join of barcode
    csvtk join -tlf 1 bam_info_cut.tsv bc_extract.tsv --left-join \
        | csvtk rename -tl -f Read,Pos,EndPos,Ref,MapQual -n read_id,start,end,chr,mapq -o read_tags.tsv

    rm bam_info.tsv bam_info_cut.tsv bc_extract.tsv
    """
}


process combine_adapter_tables {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta), path("adapters.tsv")
    output:
        tuple val(meta), path("read_config.tsv"), emit: read_config
    """
    # Concatenate tsv file keeping header from first file.
    awk 'FNR==1 && NR!=1{next;}{print}' adapters.tsv* > "read_config.tsv"
    """
}


process summarize_adapter_table {
    label "singlecell"
    cpus 1
    memory "2 GB"
    input:
        tuple val(meta), path(read_config)
    output:
        tuple val(meta), path("config_stats.json"), emit: config_stats
    """
    workflow-glue summarise_adapters \
        --read_config_tsv "${read_config}" \
        --sample_id "${meta.alias}" \
        --out "config_stats.json" \
        --threads $task.cpus
    """
}


// workflow module
workflow preprocess {
    take:
        read_chunks
        bc_longlist_dir
        ref_genome_fasta
        ref_genome_idx
        ref_genes_gtf
    main:
        // alignment pre-requisites
        call_paftools(ref_genes_gtf)
        get_chrom_sizes(ref_genome_idx)
        build_minimap_index(ref_genome_fasta)
        
        // find adapters, trim barcodes, and align
        call_adapter_scan(
            read_chunks,
            bc_longlist_dir,
            build_minimap_index.out.index,
            call_paftools.out.ref_genes_bed,
            get_chrom_sizes.out.ref_chrom_sizes)

        // TODO: this could/should be moved out of here as its not really part
        //       of the per-chunk alignment processing. For instance we aren't
        //       also combining the barcode counts here.
        combine_adapter_tables(
            call_adapter_scan.out.read_config_chunked.groupTuple())
        summarize_adapter_table(combine_adapter_tables.out.read_config)

    emit:
        bam_sort = call_adapter_scan.out.bam_sort
        read_tags = call_adapter_scan.out.read_tags
        high_qual_bc_counts = call_adapter_scan.out.barcode_counts
        config_stats = summarize_adapter_table.out.config_stats
}
