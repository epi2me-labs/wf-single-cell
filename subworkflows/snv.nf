OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")


process make_genome_dict {
    label "singlecell"
    cpus 1
    memory "4 GB"
    input:
        path "ref.fa"
        path "ref.fa.fai"
    output:
        path "ref.dict", emit: ref_genome_dict
    script:
    """
    gatk CreateSequenceDictionary --REFERENCE ref.fa
    """
}

process get_contigs {
    label "singlecell"
    cpus 2
    memory "2 GB"  
    input:
        tuple val(meta),
              path('tagged.bam'),
              path('tagged.bam.bai')
    output:
          tuple val(meta),
                path('contigs.csv')
    script:
    """
    # Get the contigs with any reads mapping
    samtools idxstats tagged.bam | awk '\$3 > 0 {print \$1}' > contigs.csv
    """
}

// UMI-deduplication and exon-splitting.
process preprocess_and_cell_snv1 {
    label "singlecell"
    cpus 1
    memory "6 GB"
    input:
        tuple val(meta),
              val(bc),
              path('tagged.bam'),
              path('tagged.bam.bai')
        path "ref.fa"
        path "ref.fa.fai"
        path "ref.dict"
    output:
        tuple val(meta),
              path("exon_split.bam"),
              path("exon_split.bam.bai"),
              val(bc),
              emit: exon_split_bam
        tuple val(meta),
              path("${bc}.vcf.gz"),
              val(bc),
              emit: vcf

    script:
    """
    # 1. Deduplicate UMIs; representative record selected that has the 
    # lowest number of mapping coordinates or highest mapping quality.
    # Use the --unique method as we have corrected UMIs earlier in the workflow
    # umi_tools dedupe does not read from stdin 
    
    # 2.Split records by exon as longshot does not like the 'N's separating them.
    # After exon splitting, rare cases of trailing Ns at the end of CIGAR (# eg. 100M99N10H)
    # appear to cause GATK SplitNCigarReads to fail
    # This may be have fixed by removing --cap-sw-mem 50m from the mm2 command

    log(){
        echo "\$(date +"%Y-%m-%d %H:%M:%S") - \$1"
    }
    
    mkdir -p tmp_dir
    
    umi_tools dedup \
        --umi-tag UB \
        --temp-dir tmp_dir \
        --cell-tag CB \
        --extract-umi-method=tag \
        --method unique \
        --per-cell -I tagged.bam -S dedup.bam
    log "UMI deduplication done"

    gatk SplitNCigarReads -R ref.fa \
        --input dedup.bam \
        --output exon_split_tmp.bam \
        --tmp-dir tmp_dir \
        --java-options "-Xmx8g" \
        --max-reads-in-memory 100000
    log "Exon splitting done"
    
    # We will be splitting by barcode later, 
    # so lets save space by removing @PG lines in the header
    samtools reheader -c 'grep -v ^@PG' exon_split_tmp.bam  > exon_split.bam
    samtools index -@ $task.cpus exon_split.bam

    rm exon_split_tmp.bam*

    # Do the first round cell SNV calling
    vcf="${bc}.vcf"
    
    log "Running longshot"
    longshot \
        --bam exon_split.bam \
        --ref ref.fa \
        --min_alt_count 1 \
        --min_cov 1 \
        --sample_id $bc \
        --out \$vcf
    bgzip \${vcf} 
    """
}

process merge_processed {
    label "singlecell"
    // This process can be really slow, and should be improved. 
    // Until then, throw lots of threads at it
    cpus params.wf.merge_threads
    memory "31 GB"
    input:
        tuple val(meta),
              path('bams/cell*.bam'),
              path('bais/cell*.bai')  
    output:
        tuple val(meta),
              path('merged.bam'),
              path("merged.bam.bai") 

    script:
    """
    samtools merge --no-PG -@ $task.cpus bams/*.bam --write-index -o merged.bam##idx##merged.bam.bai
    """
}

process split_bam_by_bc {
    label "singlecell"
    cpus params.threads
    memory "31 GB" // 13.6GB 
    input:
        tuple val(meta),
              path('tagged.bam'),
              path('tagged.bam.bai')  
    output:
        tuple val(meta),
              path('per_cell_bams/*'),
              emit: per_cell_bams 

    script:
    """
    # This is going to create lots of files!
    
    ulimit -n 20000 
    mkdir per_cell_bams

    samtools split --max-split 20000 \
            -d CB \
            --threads $task.cpus \
            -f 'per_cell_bams/%!.bam' \
            --no-PG \
            tagged.bam
    samtools index -@ $task.cpus -M per_cell_bams/*.bam
    """
}


// Initial round of bulk genotyping on the per-chr BAMs.
process snv_bulk {
    label "singlecell"
    cpus 1
    memory "16 GB"  // failiry constant ~ 8GB
    input:
        tuple val(meta),
            path('exon_split.bam'),
            path('exon_split.bam.bai'),
            val(chr)
        path "ref.fa"
        path "ref.fa.fai"
    output:
        tuple val(meta),
              path('snv_bulk.vcf.gz'),
              emit: vcf 
    script:
    """
    longshot --bam exon_split.bam \
        --ref ref.fa \
        --min_alt_count 2 \
        --min_cov 2 \
        --out snv_bulk.vcf \
        --region $chr
    bgzip snv_bulk.vcf
    """
}


process bulk_cell_merge {
    label "singlecell"
    cpus params.wf.merge_threads
    memory "31 GB"
    input:
        tuple val(meta),
              path('bulk_chunks/bulk*.vcf.gz'),
              path('per_cell_vcfs/per_cell?.vcf.gz')
    output:
        tuple val(meta),
              path('merged.vcf.gz'),
              path('merged.vcf.gz.tbi'),
              emit: merged_vcf
    script:
    merge_threads = Math.max(2, task.cpus - 2)
    """
    find -L bulk_chunks -name "*.vcf.gz" > bulk_concat.txt
    bcftools concat --file-list bulk_concat.txt --threads task.cpus --output-type z -o bulk_merged.vcf.gz
    >&2 echo "Finished merging bulk VCFs"

    tabix bulk_merged.vcf.gz
    echo bulk_merged.vcf.gz > merge_list.txt

    find -L per_cell_vcfs -name '*.vcf.gz' -exec tabix {} \\;
    find -L per_cell_vcfs -name "*.vcf.gz" >> merge_list.txt
    
    # Merge and then remove the cell columns (keep SAMPLE only)
    # now we have our merged VCF of all sites
    # This speeds up the next retrospective genotyping step
    bcftools merge --file-list merge_list.txt \
        | bcftools norm --atomize --multiallelics - \
        | bcftools +fill-tags - -- -t all  \
        | bcftools view -s SAMPLE - --threads $merge_threads --output-type z -o merged.vcf.gz
    
    tabix merged.vcf.gz
    """
}


// Second round of single-cell genotyping of the per-cell BAMs.
process snv_cell_2 {
    label "singlecell"
    cpus 1
    memory { 6.GB * task.attempt } 
    maxRetries 1
    errorStrategy { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta),
              path('bulk_candidates.vcf.gz'),
              path('bulk_candidates.vcf.gz.tbi'),
              path(cell_bam),
              path(cell_bai),
              val(barcode)
        path "ref.fa"
        path "ref.fa.fai"
    output:
        tuple val(meta),
              path('variant_debug_dir/4.0.final_genotypes.vcf.gz')
    script:
    """ 
    mkdir -p vcfs
    mkdir -p output_vcfs

    # 1. Filter candidates for each BAM where there is any BAM coverage
    # 2. Get read depth at these locations.
    # 3. Clip the max-depth at these regions to 200x    
    # 4. Genotype the depth-clipped BAMs

    log(){
        echo "\$(date +"%Y-%m-%d %H:%M:%S") - \$1"
    }

    cat ref.fa.fai | cut -f1,2 > chr.sizes
    
    vcf="${barcode}.vcf"  # This is not the the vcvf that is emitted

    log "Processing $barcode"

    # We'll use this input bam if no max read clipping could be done 
    longshot_input=$cell_bam

    # For each cell, only keep candidates that have some coverage in the respective
    # BAM. This will dramatically reduce the number of candidates to genotype.

    # Get list of candidate variants that overlap reads in the cell
    bedtools intersect -header -u -a bulk_candidates.vcf.gz -b $cell_bam | bgzip > per_cell_candidates.vcf.gz
    n_candidates=\$( zcat per_cell_candidates.vcf.gz | grep -c -v '^#' || true)
    log "\$n_candidates candidates for $barcode"

    if [ \$n_candidates -gt 0 ]; then
        tabix per_cell_candidates.vcf.gz

        # Make a BED of variant locations as targets for depth clipping
        # uniq because some sites are multi allelic and we are just interested in the position

        zgrep -v '^#' per_cell_candidates.vcf.gz | awk  'BEGIN {OFS="\t"} {print \$1, \$2, \$2}' | uniq \
            | bedtools slop -i - -g chr.sizes -b 100 \
            | bedtools merge > variant_windows.bed
    
        # If there's no coverage for this cell at the candidate sites, just output the empty VCF
        if [ -s variant_windows.bed ]; then
            # Reduce coverage to a max of 200, if the depth is 10% higher than this we don't try to clip any more
            # This needs to be done per cell as low expressed genes in some cells may get
            # wiped out otherwise
            log calculating read depth
            mkdir mosdepth
            mosdepth --no-per-base --by variant_windows.bed --threads 4 --fast-mode  mosdepth/md $cell_bam
            csvtk filter --tabs --no-header-row --filter '4>220' mosdepth/md.regions.bed.gz > high_cov.bed
            rm -rf mosdepth

            if [ ! -s high_cov.bed ]; then
                echo "No high coverage regions found for $barcode". Skipping depth clip
            else            
                n_high_cov_regions=\$(wc -l high_cov.bed | cut -f1 | xargs )
                log "\${n_high_cov_regions} high coverage regions to be depth-clipped for $barcode"

                # Get the read_ids and depth at the high coverage regions
                # There could be multiple windows for the same read_ids
                # due to supplementary reads. We will just take the first read for now
                bedtools intersect -a $cell_bam -b high_cov.bed -wb -bed \
                    | cut -f 4,16 | sort -u -k1,1 > to_clip.bed 

                workflow-glue clip_depth \
                    --bed to_clip.bed \
                    --bam_in $cell_bam \
                    --bam_out depth_clipped.bam \
                    --target_depth 200

                samtools index depth_clipped.bam
                longshot_input=depth_clipped.bam
            fi
        fi
    fi    


    # The emitted file is the penultimate VCF: 4.0.final_genotypes.vcf from the debug directory. 
    # This preserves 0/0 variants from the individual cell genotyping, 
    # as we still want to quantify representation of REF alleles in each cell.

    log "Genotyping $barcode"
    longshot \
        --potential_variants per_cell_candidates.vcf.gz \
        --bam \$longshot_input \
        --ref ref.fa \
        --min_alt_count 1 \
        --min_cov 1 \
        --sample_id "$barcode" \
        --out "\$vcf" \
        --force_overwrite \
        --variant_debug_dir variant_debug_dir
    
    # If no VCF output, because lack of candidates for example, emit the intermediate VCF
    if ! [ -f "variant_debug_dir/4.0.final_genotypes.vcf" ]; then
        mv "variant_debug_dir/1.0.potential_SNVs.vcf" "variant_debug_dir/4.0.final_genotypes.vcf"
    fi
    bgzip "variant_debug_dir/4.0.final_genotypes.vcf"
    rm -f depth_clipped.bam*
    """ 
}


process merge_cell_snv {
    label "singlecell"
    cpus params.wf.merge_threads
    memory "31 GB"
    publishDir "${params.out_dir}/${meta.alias}", 
        mode: 'copy', 
        pattern: "*{final_merged.vcf}*"
    publishDir "${params.out_dir}/${meta.alias}", 
        mode: 'copy', 
        pattern: "${meta.alias}.genotype_matrix/*"
    input:
        tuple val(meta),
              path('round_2_cell_vcfs/cell_vcf*.vcf.gz')
        path(report_variants)
    output:
        tuple val(meta),
              path("${meta.alias}.final_merged.vcf.gz"),
              path("${meta.alias}.final_merged.vcf.gz.tbi"),
              emit: final_vcf
        tuple val(meta),
              path("${meta.alias}.genotype_matrix/*"),
              emit: genotype_matrix
        tuple val(meta),
              path("top_snvs.tsv"),
              emit: top_snvs

    script:
        // Get the user variants of interest. 
        // If not an optional file, write variant IDs frome the VCF given to vars_to_report.tsv
        opt_report_variants = report_variants.fileName.name == "OPTIONAL_FILE" ?  "" : "--report_vars vars_to_report.tsv"
        merge_threads = Math.max(2, task.cpus - 2)
    """
    # Merge the single cell VCFs
    # Use a sequential merging strategy to avoid memory issues and speed things up
    # https://bioinformatics.stackexchange.com/a/20900

    log(){
        echo "\$(date +"%Y-%m-%d %H:%M:%S") - \$1"
    }

    # merge options
    # --atomize: split MNVs into consecutive SNVs. 
    # --multiallelics - : split multiallelic sites into biallelic records
    # fill-tags - -- -t all: Fills in missing tags

    chunk_size=200
    mkdir merge_list_chunks
    mkdir merged_chunks

    log "Indexing VCFs"
    ulimit -n 20000
    find -L round_2_cell_vcfs -name '*.vcf.gz' -exec tabix {} \\;
    log "Getting VCF file list"
    find -L round_2_cell_vcfs -name "*.vcf.gz" -exec echo {} \\; >> merge_list.txt

    if [ \$(wc -l  < merge_list.txt) -gt \$chunk_size ]; then
        log "splitting file list" 
        split --lines \$chunk_size --additional-suffix '.txt' merge_list.txt merge_list_chunks/
    
        # Merge the subsets
        for sub in merge_list_chunks/*; do
            name=\$(basename \$sub)
            echo "Merging subset \${name}"
            out_vcf="merged_chunks/\${name}.vcf.gz"
            bcftools merge --file-list "\$sub" -o "\$out_vcf" --threads $merge_threads
            tabix "\$out_vcf"
        done
        find merged_chunks -name "*.vcf.gz" > final_merge_list.txt
    else
        mv merge_list.txt final_merge_list.txt
    fi

    # Final merge of the merged subsets    
    bcftools merge --file-list final_merge_list.txt \
        | bcftools norm --atomize --multiallelics - \
        | bcftools +fill-tags - -- -t all \
        | bgzip -@ $merge_threads > "${meta.alias}.final_merged.vcf.gz"
    tabix "${meta.alias}.final_merged.vcf.gz" 

    # Make a sparse MTX format snv x cell matrix where:
    # hom ref: 0
    # het: 1
    # hom alt: 2

    # If we have a list of requested variants locations, extract these from the matrix
    if ! [ -z "$opt_report_variants" ]; then
        # Get variants of interest for report. Limit to 50
        log "Getting variants of interest"
        bcftools query --format '%CHROM\\_%POS' $report_variants  | head -n 50 > vars_to_report.tsv
    fi

    # Get sparse matrix output to stdout
    mkdir "${meta.alias}.genotype_matrix"
    workflow-glue variant_mex \
        "${meta.alias}.final_merged.vcf.gz" "${meta.alias}.genotype_matrix" $opt_report_variants \
        | bgzip -@4 >  "matrix.mtx.gz"

    # cat the header and the matrix as we don't know header dims upfront
    bgzip header.txt
    cat header.txt.gz matrix.mtx.gz > "${meta.alias}.genotype_matrix/matrix.mtx.gz"
    rm matrix.mtx.gz
    """
}


workflow longshot {
    take: 
        tagged_bam
        ref_genome_fasta
        ref_genome_idx
    main:

        if (params.report_variants){
            report_variants = file("${params.report_variants}", checkIfExists: true)
        }else {
            report_variants = OPTIONAL_FILE
        }

        make_genome_dict(ref_genome_fasta, ref_genome_idx)
        
        get_contigs(tagged_bam)
        
        split_bam_by_bc(tagged_bam)

        split_bam_by_bc.out.per_cell_bams
            .flatMap{meta, files -> 
                // get [meta, barcode, bam or bai]
                files.collect { file -> [meta, file.getName().split( /\./ )[0] , file] } }
            .branch {
                meta, bc, file ->
                bam : file.name.endsWith('.bam')
                bai: file.name.endsWith('.bai')
            }.set { cell_bams }
        
        // [meta, barcode, bam, bai]
        cells = cell_bams.bam.join(cell_bams.bai, by: [0, 1])

        // First round of SNV calling per cell
        // This will identify candidate SNVs that are well represented in individual cells
        // that might be lost amongst the bulk data
        preprocess_and_cell_snv1(cells,
            ref_genome_fasta, ref_genome_idx, make_genome_dict.out)
        
        merge_processed(
            preprocess_and_cell_snv1.out.exon_split_bam
                .groupTuple()
                .map {meta, bams, bais, _bcs -> [meta, bams, bais]})
        
        // Identify candidate SNVs in bulk, this will identify candidates that
        // would not be identified from single-cell due to low coverge per cell 
        // Split bulk Bam by chromosome to speed things up
        snv_bulk(
            merge_processed.out.combine(
                get_contigs.out
                    .splitCsv(header:['chr'])
                    .map {meta, row -> [meta, row.chr]},
                by:0), 
            ref_genome_fasta, ref_genome_idx)
        
        // Merge the round 1 single-cell and bulk variants. 
        bulk_cell_merge(
            snv_bulk.out.vcf.groupTuple()
                .join(preprocess_and_cell_snv1.out.vcf.groupTuple())
                .map {meta, bulk_vcfs, cell_vcfs, _barcodes -> [meta, bulk_vcfs, cell_vcfs.flatten()]})


        // Call variants on the per-cell BAM chunks using the candidate variants 
        // from the merged bulk+cell VCFs 
        snv_cell_2(
            bulk_cell_merge.out.merged_vcf
                // [meta, merged_vcf, .tbi, bam, bai, BC]
                .combine(preprocess_and_cell_snv1.out.exon_split_bam.groupTuple().transpose(), by: 0),
            ref_genome_fasta, ref_genome_idx)

        merge_cell_snv(snv_cell_2.out.groupTuple(), report_variants)
    
    emit: top_snvs = merge_cell_snv.out.top_snvs
        
}