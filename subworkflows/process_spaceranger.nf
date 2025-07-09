include { call_paftools; build_minimap_index} from '../modules/local/common'

// Move to common
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

process get_bam_chrs {
    label "singlecell"
    cpus 1
    memory "4 GB"
    input:
        tuple val(meta),
            path("sr.bam"),
            path("sr.bam.bai")
    output:
        tuple val (meta),
              path("chr_list.txt")
    script:
    """
    samtools view -H sr.bam | grep '^@SQ' | cut -f2 | sed 's/SN://' \
        > chr_list.txt    
    echo '*' >> chr_list.txt
    """
}


process process_long_reads {
    label "singlecell"
    cpus params.threads
    memory 17.GB  
    input:
        tuple val(meta), 
              path(chunk, stageAs: 'chunk.fq.gz')
        path('genome_index.mmi')
        path('ref_genes.bed')

    output:
        tuple val(meta),
              path("read_tags.tsv"),
              emit: mapping_tags
        tuple val(meta),
              path("sorted.bam"), path("sorted.bam.bai"),
              emit: bam_sort
        tuple val(meta),
              path("bamstats.tsv"),
              emit: bam_stats
    script:
    def mm2_threads = task.cpus - 1
    """
    # -uf, these are stranded in mRNA sense ?

    # Extract fastq with relevant tags from spaceranger
    # Map reads to the reference genome, carry over the tags
    
    # Note: the minimap command is the same as that in preprocess.nf:call_adapter_scan
    # Can remove this dupoication

    minimap2 -y -ax splice -uf --MD \
        -t $mm2_threads -K 10M \
        --junc-bed ref_genes.bed  \
        --cap-kalloc 100m \
        genome_index.mmi chunk.fq.gz \
    | samtools view -uh --no-PG - \
    | tee >(seqkit bam -s  2> bamstats.tsv ) \
    | tee >(samtools view - -d SA \
        | awk 'BEGIN{OFS="\t"; print "read_id", "SA"} {print \$1,"True"}' > SA_tags.tsv ) \
    | samtools view -uh -F 256 - \
    | tee >(samtools sort --write-index -o "sorted.bam"##idx##"sorted.bam.bai" --no-PG  -) \
    | tee >(seqkit bam -F - 2> bam_info.tsv)

    csvtk cut -tlf Read,Pos,EndPos,Ref,MapQual bam_info.tsv \
    | csvtk rename -tl -f Read,Pos,EndPos,Ref,MapQual -n read_id,start,end,chr,mapq -o read_tags_interim.tsv

    # Merge the SA column with the read tags on read_id
    if [ \$(wc -l < SA_tags.tsv) -eq 1 ]; then
        echo "No SA tags found"
        # Add an empty SA column
        csvtk mutate2 -t -n 'SA' -e " '' " read_tags_interim.tsv > read_tags.tsv
    else
        csvtk -t uniq SA_tags.tsv | csvtk join -t --left-join --fields read_id read_tags_interim.tsv - > read_tags.tsv
    fi
    # rm bam_info.tsv bam_info_cut.tsv bc_extract.tsv read_tags_interim.tsv
    """    
}


process parse_sr_bam {
    // Extract barcode counts from teh spaceranger BAM
    label "wf_common"
    cpus params.threads
    memory "16 GB"
    input:
        tuple val(meta),
              path('spaceranger.bam'),
              path('spaceranger.bam.bai'),
              val(chr)
    output:
        tuple val(meta),
              path("spaceranger_tags.tsv"),
              emit: spaceranger_tags
        tuple val(meta),
              path("barcode_counts.tsv"),
              emit: barcode_counts
    script:
    """
    # If we want qual (UY, UB) these must be quoted
    workflow-glue tags_from_bam \
        spaceranger.bam \
        spaceranger_tags.tsv \
        barcode_counts.tsv \
        --tags CR CB CY UR UB UY \
        --chrom "${chr}" 
    """
}

// Aggregate the per-chunk cell counts into a per sample/chr counts file.
process aggregate_barcode_counts {
    label "wf_common"
    cpus params.threads
    memory "32 GB"
    input:
        tuple val(meta),
              path('barcode_counts/bc_counts_*.tsv')
    output:   
        path("agg_barcode_counts.tsv"),
              emit: barcode_counts
        tuple val(meta),
              path("barcode_list.tsv"),
              emit: barcode_list
    script:
    """
    csvtk concat -t barcode_counts/*.tsv \
    | csvtk -t sort -k barcode - \
    | awk -F'\t' '
        NR==1 { print; next }
        {
        if (\$1 == prev) {
            sum += \$2
        } else {
            if (NR > 2) print prev, sum
            prev = \$1
            sum = \$2
        }
        }
        END {
        if (NR > 1) print prev, sum
        }' OFS='\t' \
    | csvtk -t mutate2 --name 'sample' --expression '"${meta.alias}"' \
    | csvtk -t sort -k count:n  > agg_barcode_counts.tsv


    # Get whitelist. This is just the list of barcodes. This is identical to the list
    # in cell_counts.tsv, but this is need to be compatiable with other 10x kits
    # where filtering of cells occurs
    tail +2 agg_barcode_counts.tsv | cut -f 1 > barcode_list.tsv
    """
}


// Put in common?
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


process combine_mapping_and_demux_tags {
    // Join the demux and mapping tags and split by chromosome
    label "wf_common"
    cpus 6
    memory "30 GB" 
    input:
        tuple val(meta),
              path('demux_tags/demux_??.tsv'),
              path('mapping_tags/mapping_??.tsv')
    output:
        tuple val(meta),
              path("chr_tags/*"),
              emit: chr_tags
    script:
    def mem_mb = task.memory.toMega().toDouble()
    def buffer_mb = (Math.max(4096.0, mem_mb - 4096.0) / 2).toInteger()
    def buffer_size = "${buffer_mb}M"
    def threads_per_sort = Math.max(1, (int)(task.cpus / 2))
    def tmpdir = "tmpdir"
    """
    mkdir ${tmpdir}
    export TMPDIR=${tmpdir}
    
    # Function to concatenate + sort multiple TSVs by first column
    function concat_and_sort() {
        local input_dir="\$1"
        local output_file="\$2"
    
        # Get header from the first file
        head -n 1 "\$(find -L \"\$input_dir\" -type f -name '*.tsv' | head -n 1)" > "\$output_file"
    
        # Concatenate bodies and pipe directly to sort
        find -L "\$input_dir" -type f -name '*.tsv' | while read -r file; do
            tail -n +2 "\$file"
        done | sort \\
            --buffer-size=${buffer_size} \\
            --temporary-directory=${tmpdir} \\
            --parallel=${threads_per_sort} \\
            --field-separator=\$'\\t' \\
            --key=1,1 \\
            >> "\$output_file"
    }
    
    # Combine mapping and demux tags
    concat_and_sort demux_tags demux.tsv &
    concat_and_sort mapping_tags map.tsv &
    wait
    rm -rf ${tmpdir}/*
    
    # Join tags and split by chromosome
    # TODO: just do the splitting in the Python script?
    mkdir -p chr_tags
    workflow-glue join_tags map.tsv demux.tsv | \\
        csvtk split \\
            --num-cpus ${task.cpus} \\
            --tabs \\
            --fields chr \\
            --out-file chr_tags
    rm -f map.tsv demux.tsv

    # Clean up filenames: remove 'stdin-' prefix from each output
    shopt -s nullglob
    for f in chr_tags/stdin-*.tsv; do
        mv "\$f" "chr_tags/\${f##*/stdin-}"
    done
    
    # Find the column number for 'CB'
    A_FILE=\$(find chr_tags -type f -name '*.tsv' | head -n 1)
    CB_COL=\$(head -n 1 \${A_FILE} | awk -F'\\t' '{for (i=1; i<=NF; i++) if (\$i=="CB") print i; exit}')
    
    # Sort each per-chromosome file by CB
    for file in chr_tags/*.tsv; do
      head -n 1 "\$file" > "\${file}.sorted"
      tail -n +2 "\$file" | sort \\
        --buffer-size=${buffer_size} \\
        --temporary-directory=${tmpdir} \\
        --parallel=${task.cpus} \\
        --field-separator=\$'\\t' \\
        --key=\${CB_COL},\${CB_COL} \\
        >> "\${file}.sorted"
      mv "\${file}.sorted" "\$file"
    done
    
    rm -rf ${tmpdir}
    """
}


workflow spaceranger {
    take:
        long_reads
        spaceranger_bam
        ref_genome_fa
        ref_genes_gtf

    main:
        ref_genome_mmi = build_minimap_index(ref_genome_fa)
        ref_genome_bed = call_paftools(ref_genes_gtf)

        get_bam_chrs(spaceranger_bam)
        spaceranger_chrs = get_bam_chrs.out
            .splitCsv(header:false)

        chr_gtfs = split_gtf_by_chroms(ref_genes_gtf)
            .flatten()
            .map {fname -> tuple(fname.baseName, fname)} // [chr, gtf]

        adapter_stats = long_reads.groupTuple()
                .map {meta, _fastq -> [ meta, file("${meta.adapter_stats}", checkIfExists: true)]}
        
        process_long_reads(
            long_reads, ref_genome_mmi, ref_genome_bed)
        
        merged_long_read_bam = merge_bams(process_long_reads.out.bam_sort.groupTuple())
        
        parse_sr_bam(spaceranger_bam
            .combine(
                spaceranger_chrs.map { meta, chrs -> [meta, chrs[0]] }, by: 0))
        
        // Combine BC counts by sample
        aggregate_barcode_counts(
            parse_sr_bam.out.barcode_counts.groupTuple())

        // Merge demux tags with features and cat by chromosome
        combine_mapping_and_demux_tags(
            parse_sr_bam.out.spaceranger_tags.groupTuple()
            .join(process_long_reads.out.mapping_tags.groupTuple()))

        // Extract chromosome name from the tags filenames, vom
        chr_tags = combine_mapping_and_demux_tags.out.chr_tags
        .transpose()
            .map {meta, tags -> 
                def chr = tags.baseName
                [meta, chr, tags]
            }

    emit:
        adapter_summary = adapter_stats
        chr_tags = chr_tags
        merged_bam = merged_long_read_bam
        bam_stats = process_long_reads.out.bam_stats
        barcode_list = aggregate_barcode_counts.out.barcode_list
        barcode_counts = aggregate_barcode_counts.out.barcode_counts
            .collectFile(keepHeader: true)
        
}
