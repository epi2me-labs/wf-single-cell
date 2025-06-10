include { call_paftools; build_minimap_index} from '../modules/local/common'

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
    output:
        tuple val(meta), path("adapters.json"), emit: adapter_summary
        tuple val(meta), path("read_tags.tsv"), emit: read_tags
        tuple val(meta), path("high_quality_bc_counts.tsv"), emit: barcode_counts
        tuple val(meta), path("sorted.bam"), path("sorted.bam.bai"), emit: bam_sort
        tuple val(meta), path("bamstats.tsv"), emit: bam_stats
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
        --summary "adapters.json" \
        ${fl} \
    | workflow-glue extract_barcode \
        - \
        bc_longlist_dir/${meta['bc_long_list']} \
        --kit ${meta["kit_name"]} \
        --adapter1_suff_length $params.barcode_adapter1_suff_length \
        --min_barcode_qv $params.barcode_min_quality \
        --barcode_length ${meta['barcode_length']} \
        --umi_length ${meta['umi_length']} \
        --output_read_tags "bc_extract.tsv" \
        --output_barcode_counts "high_quality_bc_counts.tsv" \
    | minimap2 -ax splice -uf --MD \
        -t $mm2_threads -K 10M \
        --junc-bed ref_genes.bed  \
        --cap-kalloc 100m \
        genome_index.mmi - \
    | samtools view -uh --no-PG - \
    | tee >(seqkit bam -s  2> bamstats.tsv ) \
    | tee >(samtools view - -d SA \
        | awk 'BEGIN{OFS="\t"; print "read_id", "SA"} {print \$1,"True"}' > SA_tags.tsv ) \
    | samtools view -uh -F 256 - \
    | tee >(samtools sort --write-index -o "sorted.bam"##idx##"sorted.bam.bai" --no-PG  -) \
    | seqkit bam -F - 2> bam_info.tsv

    # TODO: improve this with pipes?
    csvtk cut -tlf Read,Pos,EndPos,Ref,MapQual bam_info.tsv > bam_info_cut.tsv
    # Left join of barcode
    csvtk join -tlf 1 bam_info_cut.tsv bc_extract.tsv --left-join \
        | csvtk rename -tl -f Read,Pos,EndPos,Ref,MapQual -n read_id,start,end,chr,mapq -o read_tags_interim.tsv

    # Merge the SA column with the read tags on read_id
    if [ \$(wc -l < SA_tags.tsv) -eq 1 ]; then
        echo "No SA tags found"
        # Add an empty SA column
        csvtk mutate2 -t -n 'SA' -e " '' " read_tags_interim.tsv > read_tags.tsv
    else
        csvtk -t uniq SA_tags.tsv | csvtk join -t --left-join --fields read_id read_tags_interim.tsv - > read_tags.tsv
    fi
    rm bam_info.tsv bam_info_cut.tsv bc_extract.tsv read_tags_interim.tsv
    """
}


process merge_bams {
    // Combine all BAMs derived from the initial chunking into per sample files
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
            path('bams/*aln.bam'),
            path('bams/*aln.bam.bai')
    output:
        tuple val(meta),
              path("merged.sorted.bam"),
              path("merged.sorted.bam.bai"),
              emit: merged_bam
    script:
    """
    samtools merge -@ ${task.cpus -1} --write-index -o "merged.sorted.bam##idx##merged.sorted.bam.bai" bams/*.bam
    """
}


// workflow module
workflow preprocess {
    take:
        read_chunks
        bc_longlist_dir
        ref_genome_fasta
        ref_genes_gtf
    main:
        // alignment pre-requisites
        index_mmi = build_minimap_index(ref_genome_fasta)
        ref_genes_bed = call_paftools(ref_genes_gtf)


        // find adapters, trim barcodes, and align
        call_adapter_scan(
            read_chunks,
            bc_longlist_dir,
            build_minimap_index.out.index,
            ref_genes_bed)
        
        merged_bam = merge_bams(call_adapter_scan.out.bam_sort.groupTuple())

    emit:
        merged_bam = merged_bam
        bam_stats = call_adapter_scan.out.bam_stats
        read_tags = call_adapter_scan.out.read_tags
        high_qual_bc_counts = call_adapter_scan.out.barcode_counts
        adapter_summary = call_adapter_scan.out.adapter_summary
}
