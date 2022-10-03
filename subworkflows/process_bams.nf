import java.util.ArrayList;

process get_kit_info {
    label "singlecell"
    cpus 1
    input:
        path kit_config
        path sc_sample_sheet
              
    output:
        path 'sample_kit_info.csv', emit: kit_info
    
    script:
    """
    #!/usr/bin/env python
    import pandas as pd

    sample_df = pd.read_csv("${sc_sample_sheet}", sep=",", comment="#").set_index("sample_id", drop=True)

    kit_df = pd.read_csv("$kit_config", sep=",", comment="#")


    records = []
    for sample_id, row in sample_df.iterrows():
        kit_name = sample_df.loc[sample_id, "kit_name"]
        kit_version = sample_df.loc[sample_id, "kit_version"]

        # Get barcode length based on the kit_name and kit_version specified for this sample_id.
        rows = (kit_df["kit_name"] == kit_name) & (kit_df["kit_version"] == kit_version)
        barcode_length = kit_df.loc[rows, "barcode_length"].values[0]

        # Get the appropriate cell barcode longlist based on the kit_name specified for this sample_id.
        if kit_name == "3prime":
            if kit_version == "v2":
                long_list = "737K-august-2016.txt.gz"
            elif kit_version == "v3":
                long_list = "3M-february-2018.txt.gz"
            else:
                raise Exception("Encountered an unexpected kit version for 3prime kit (v2 or v3)")
        elif kit_name == "5prime":
            long_list = "737K-august-2016.txt.gz"
        elif kit_name == "multiome":
            long_list = "737K-arc-v1.txt.gz"
        else:
            raise Exception("Encountered an unexpected kit_name in samples.csv")

        # Get UMI length based on the kit_name and kit_version specified for this sample_id.
        umi_length = kit_df.loc[rows, "umi_length"].values[0]
        records.append([sample_id, kit_name, kit_version, barcode_length, umi_length, long_list])
    df_out = pd.DataFrame.from_records(records)
    df_out.columns = ["sample_id", "kit_name", "kit_version", "barcode_length", "umi_length", "long_list"]
    df_out.to_csv('sample_kit_info.csv', index=False)
    """
}


process extract_barcodes{
    /*
    Build minimap index from reference genome
    */
    label "singlecell"
    cpus params.max_threads
    input:
        tuple val(sample_id), 
              val(kit_name),
              val(kit_version),
              val(barcode_length),
              val(umi_length),
              val(bc_long_list),
              path(bam_sort),
              path(bam_sort_idx)
        path bc_longlist_dir

    output:
        tuple val(sample_id), 
              path("*.bc_extract.sorted.bam"), 
              path("*.bc_extract.sorted.bam.bai"), 
              emit: bam_bc_uncorr
        tuple val(sample_id), 
              path("*.tsv"), 
              emit: barcode_counts
    """
    extract_barcode.py \
    $bam_sort ${bc_longlist_dir}/${bc_long_list}\
    -t $task.cpus \
    --kit $kit_name \
    --adapter1_suff_length $params.barcode_adapter1_suff_length \
    --min_barcode_qv $params.barcode_min_quality \
    --barcode_length $barcode_length \
    --umi_length $umi_length \
    --output_bam "tmp.bam" \
    --output_barcodes "${sample_id}.uncorrected_bc_counts.tsv";

    # Cleanup headers
    samtools reheader --no-PG -c 'grep -v ^@PG' tmp.bam > "${sample_id}.bc_extract.sorted.bam"
    samtools index "${sample_id}.bc_extract.sorted.bam"
    rm tmp.bam
    """
}

process generate_whitelist{
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              path(counts),
              val(expected_cells)
    output:
        tuple val(sample_id), path("*whitelist.tsv"), emit: whitelist
        tuple val(sample_id), path("*kneeplot.png"), emit: kneeplot
    """
    knee_plot.py \
        $counts \
        --exp_cells $expected_cells \
        --output_whitelist "${sample_id}.whitelist.tsv" \
        --output_plot "${sample_id}.kneeplot.png"
    """
}

process split_bam_by_chroms{
    label "singlecell"
    cpus params.max_threads
    input:
        tuple val(sample_id), path(bam), path(bai)
    output:
        tuple val(sample_id), 
              path("splits/*.bam"),
              path("splits/*.bai"), emit: bam
    """
    split_bam_by_chroms.py -t ${task.cpus} --output_dir splits $bam
    """
}

process split_fasta_by_chroms {
    label "singlecell"
    cpus params.max_threads
    input:
        path(ref_genome)
    output:
        path("fasta/*")
    """
    samtools faidx $ref_genome
    seqkit split $ref_genome --by-id --two-pass -O fasta
    cd fasta
    for file in *; do
        newname=\${file#*part_};
        mv \$file \$newname;
    done
    """
}


process assign_barcodes{
    label "singlecell"
    cpus 1
    input:
         tuple val(sample_id), 
               val(kit_name),
               val(kit_version),
               val(barcode_length),
               val(umi_length),
               val(bc_long_list),
               path(whitelist),
               val(_), // Redundant sample_id: remove
               val(chr), 
               path(bam),
               path(bai)
    output:
        tuple val(sample_id), 
            val(chr),
            path("*.bc_assign.bam"),
            path("*.bc_assign.bam.bai"),
            emit: chrom_bam_bc_bai
        tuple val(sample_id),
              val(chr),
              path("*.bc_assign_counts.tsv"), 
              emit: chrom_assigned_barcode_counts
    """
    assign_barcodes.py -t 1 \
        --output_bam tmp.bam \
        --output_counts ${sample_id}_${chr}.bc_assign_counts.tsv \
        --max_ed $params.barcode_max_ed \
        --min_ed_diff $params.barcode_min_ed_diff \
        --kit $kit_name \
        --adapter1_suff_length $params.barcode_adapter1_suff_length \
        --barcode_length $barcode_length \
        --umi_length $umi_length \
        $bam $whitelist
    
    samtools reheader --no-PG -c 'grep -v ^@PG' tmp.bam \
        > "${sample_id}_${chr}.bc_assign.bam";
    samtools index "${sample_id}_${chr}.bc_assign.bam"
    rm tmp.bam

"""
}


process bam_to_bed {
    label "singlecell"
    cpus 1
    input:
        tuple val(chr),
              val(sample_id),
              path(bam),
              path(bai) 
    output:
        tuple val(sample_id), 
              val(chr),
              path('*bc_assign.bed'), 
              emit: chrom_bed_bc
    """
    bedtools bamtobed -i $bam > "${sample_id}_bc_assign.bed"
    """
}

process split_gtf_by_chroms {
    label "singlecell"
    cpus 1
    input:
        path(gtf)
    output:
        path("*"), emit: chrom_gtf
    """
    gawk '/^[^#]/ {print>\$1".gtf"}' $gtf 
    """
}   

process assign_genes {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path(chrom_bed_bc),
              path(chrom_gtf)
    output:
        tuple val(sample_id),
              val(chr),
              path("*.read.gene_assigns.tsv"),
              emit: chrom_tsv_gene_assigns
    """
    assign_genes.py \
    --output ${sample_id}_${chr}.read.gene_assigns.tsv \
    $chrom_bed_bc $chrom_gtf
    """
}

process add_gene_tags_to_bam {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              val(chr),
              path(chrom_bam_bc),
              path(chrom_bam_bc_bai),
              path(chrom_tsv_gene_assigns)
    output:
        tuple val(sample_id),
              val(chr),
              path("*bc_assign.gene.bam"), 
              path("*bc_assign.gene.bam.bai"),
              emit: chrom_bam_bc_bai
    """
    add_gene_tags.py \
        --output tmp.bam \
        $chrom_bam_bc $chrom_tsv_gene_assigns
    
    samtools reheader --no-PG -c 'grep -v ^@PG' tmp.bam \
        > ${sample_id}_${chr}_bc_assign.gene.bam;
    samtools index ${sample_id}_${chr}_bc_assign.gene.bam
    """
}

process cluster_umis {
    label "singlecell"
    cpus params.umi_cluster_max_threads
    input:
        tuple val(sample_id),
              val(chr),
              path(bam),
              path(bai)
    output:
         tuple val(sample_id),
              path("*.tagged.bam"),
              path("*.tagged.bam.bai"),
              emit: bam_bc_bai
    """
    cluster_umis.py \
    $bam \
    --threads $task.cpus \
    --ref_interval $params.umi_genomic_interval \
    --cell_gene_max_reads $params.umi_cell_gene_max_reads \
    --output tmp.bam 

    samtools reheader --no-PG -c 'grep -v ^@PG' tmp.bam \
        >  ${sample_id}_${chr}.tagged.bam;
    samtools index ${sample_id}_${chr}.tagged.bam
    rm tmp.bam
    """
}


process combine_chrom_bams {
    // Merge all chromosome bams by sample_id
    label "singlecell"
    cpus params.max_threads
    input:
        tuple val(sample_id), 
              path(bams),
              path(bais)
    output:
        tuple val(sample_id), 
              path("*tagged.sorted.bam"), 
              path("*tagged.sorted.bam.bai"),
              emit: bam_fully_tagged
    """
    samtools merge -@ ${task.cpus} -o "${sample_id}.tagged.sorted.bam" $bams; 
    samtools index -@ ${task.cpus} "${sample_id}.tagged.sorted.bam";
    """
}

process preprocess_bams_for_stringtie {
    // A termporary workaround for the reads being in the reverse orientation
    // Also split bam by barcode, keeping only the largest read per umi
    // TODO: Create a consensus instead
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              path(bam), 
              path(bai)
    output:
        tuple val(sample_id),
              path('*.bam'),
              emit: flipped_bc_bams
    """
    process_bam_for_stringtie.py $bam
    """
}

process split_bam_by_barcode {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              val(chr),
              path(bam),
    output:
        tuple val(sample_id), 
              path("*.sam"), 
              emit: sams
    """
    # Split by corrected barcode (CB:Z). 
    samtools view --no-header ${bam} | gawk 'match(\$0, /CB:Z:(\\w+)/, m) {print \$0> m[1]".sam"}'
    # Prepend header
    header=\$(samtools view -H ${bam})
    for f in *.sam; do
        echo "\$header" > tmp;
        cat \$f >> tmp;
        cp tmp \$f;
    done;
    """
}

process stringtie {
    label "singlecell"
    cpus 1
    input:
        path(ref_gtf)
        tuple val(sample_id), 
              path(bam)
    output:
        tuple val(sample_id), 
              path("barcode.tmap"), 
              emit: tmap
    """
    samtools sort $bam -o sorted.bam 
    stringtie -L -v -p ${task.cpus}  -G ${ref_gtf} -l test \
        -o stringtie.gff \
        sorted.bam
  
    barcode=\$(echo $bam | sed 's/.bam//')
    gffcompare -o gffcompare -r $ref_gtf stringtie.gff
    # Add barcode column to tmap file
    sed "1s/\$/\tbarcode/; 1 ! s/\$/\t\${barcode}/" gffcompare.stringtie.gff.tmap \
        > barcode.tmap
    """
}

process build_transcript_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id), 
              path(gffcompare_tmap),
              path(gene_matrix)
    output:
        tuple val(sample_id), 
              path("*matrix.tsv"), 
              emit: transcript_matrix
    """
    # remove the headers that were mixed in.
    cat $gffcompare_tmap | gawk '/^[^ref_gene_id]/' > tmap_cleaned

    transcript_expression_matrix.py \
        --min_genes $params.matrix_min_genes \
        --min_cells $params.matrix_min_cells \
        --max_mito $params.matrix_max_mito \
        --norm_count $params.matrix_norm_count \
        --gene_matrix $gene_matrix \
        --output ${sample_id}_processed_transcript_matrix.tsv \
        tmap_cleaned
    """
}


process count_cell_gene_umi_reads {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(bam),
              path(bai)
    output:
        tuple val(sample_id),
              path("${sample_id}_cell_umi_gene.tsv"),
              emit: cell_umi_gene_tsv
    """
    cell_umi_gene_table.py \
        --output  ${sample_id}_cell_umi_gene.tsv $bam
    """
}

process umi_gene_saturation {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(cell_umi_gene_tsv)
    output:
        tuple val(sample_id),
              path("*saturation_curves.png")
    """
    calc_saturation.py \
        --output ${sample_id}.saturation_curves.png \
        $cell_umi_gene_tsv
    """
}

process construct_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(bam),
              path(bai)
    output:
        tuple val(sample_id), 
              path("*gene_expression.counts.tsv"),
              emit: matrix_counts_tsv
    """
    gene_expression.py \
        --output ${sample_id}.gene_expression.counts.tsv $bam

    """
}

process process_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_counts_tsv)
    output:
        tuple val(sample_id), 
              path("*gene_expression.processed.tsv"),
              emit: matrix_processed_tsv
        tuple val(sample_id),
              path("*gene_expression.mito.tsv"),
              emit: matrix_mito_tsv
    """
    process_matrix.py \
    --min_genes $params.matrix_min_genes \
    --min_cells $params.matrix_min_cells \
    --max_mito $params.matrix_max_mito \
    --norm_count $params.matrix_norm_count \
    --output ${sample_id}.gene_expression.processed.tsv \
    --mito_output ${sample_id}.gene_expression.mito.tsv \
    $matrix_counts_tsv

    """
}

process umap_reduce_expression_matrix {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_processed_tsv)
    output:
         tuple val(sample_id),
              path("*gene_expression.umap.tsv"),
              emit: matrix_umap_tsv
    """
    umap_reduce.py \
        --output ${sample_id}.gene_expression.umap.tsv \
        $matrix_processed_tsv
    """
}


process umap_plot_total_umis {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_processed_tsv)
    output:
          tuple val(sample_id),
              path("*umap.total.png"), 
              emit: matrix_umap_plot_total

    """
   plot_umap.py \
        --output ${sample_id}.umap.total.png \
        $matrix_umap_tsv $matrix_processed_tsv

    """
}

process umap_plot_genes {
    // TODO: make a channle of input genes for thes process
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_processed_tsv),
              val(gene)
    output:
        tuple val(sample_id),
              path("*umap.gene.${gene}.png"), 
              emit: matrix_umap_plot_gene

    """
    plot_umap.py \
        --gene $gene \
        --output ${sample_id}.umap.gene.${gene}.png \
        $matrix_umap_tsv $matrix_processed_tsv

    """
}

process umap_plot_mito_genes {
    label "singlecell"
    cpus 1
    input:
        tuple val(sample_id),
              path(matrix_umap_tsv),
              path(matrix_mito_tsv)
    output:
        tuple val(sample_id),
              path("*umap.mitochondrial.png"), 
              emit: matrix_umap_plot_mito
    """
    plot_umap.py \
        --mito_genes \
        --output ${sample_id}.umap.mitochondrial.png \
        $matrix_umap_tsv \
        $matrix_mito_tsv
    """
}
    
workflow process_bams {
    take:
        bam
        bam_idx
        sc_sample_sheet
        kit_config
        gtf
        umap_genes
        bc_longlist_dir
        sample_kits
        ref_genome_fasta
   
    main:

        get_kit_info(
            kit_config,
            sc_sample_sheet)

        get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(bam).join(bam_idx)

        extract_barcodes(
            get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(bam).join(bam_idx),
            bc_longlist_dir)
        
        split_bam_by_chroms(
            extract_barcodes.out.bam_bc_uncorr
        )

        split_fasta_by_chroms(
            ref_genome_fasta
        )

        generate_whitelist(
            extract_barcodes.out.barcode_counts
            .join(sample_kits)
            .map{it -> tuple(it[0], it[1], it[4] )}
        )
        
        // Extract chr from filename and add to tuple to give: 
        // [sample_id, chr, bam, bai]
        bam_bai_chromes = split_bam_by_chroms.out.bam.map({it ->
            sbi = []
            if (it[1].getClass() == java.util.ArrayList){
                // Multiple chroms:
                // [sample_id, [bam1, bam2], [bai1, bai2]]
                for (i=0; i<it[1].size(); i++) {
                    chr = it[1][i].name.strip('sorted.bam')
                    sbi.add(tuple(it[0], chr, it[1][i], it[2][i]))
                }
            }
            else{
                // Only single chrom so we have:
                // [sample_id, bam, bai]
                chr = it[1].name.strip('sorted.bam')
                sbi.add(tuple(it[0], chr, it[1], it[2])) 
            }
            return sbi
            }).flatMap(it-> it)
        
        // merge kit info to bams
        chr_bam_kit = get_kit_info.out.kit_info
            .splitCsv(header:false, skip:1)
            .join(generate_whitelist.out.whitelist)
            .cross(bam_bai_chromes).map({it -> it.flatten()})

        assign_barcodes(chr_bam_kit)

        bam_to_bed(assign_barcodes.out.chrom_bam_bc_bai)

        chr_gtf = split_gtf_by_chroms(gtf)
            .flatten()
            .map {file -> 
                // create [chr, gtf]
                tuple(file.baseName, file)}
        
        // combine all chr bams with chr gtfs
        chr_beds_gtf = chr_gtf.cross(
            bam_to_bed.out.chrom_bed_bc)
            .map({it ->
            //  rejig the tuple to [sample_id, chr, bed, gtf]
             tuple(it[1][1], it[0][0], it[1][2], it[0][1])})
        

        // Make tuples: sample_id, chr, bam, bai, ref_gtf, ref_genome 
        split_fasta_by_chroms.out
        chr_gtf_bam = chr_gtf.join(
            split_fasta_by_chroms.out.flatten()
            .map(it -> tuple(it.baseName, it))
            )
            .cross(
                 assign_barcodes.out.chrom_bam_bc_bai
                 .map(it -> [it[1], it[0], it[2], it[3]])
            )
            // Rejig for next process
            .map(it -> [it[1][1], it[0][0], it[1][2], it[1][3], it[0][1], it[0][2]])

        assign_genes(chr_beds_gtf)

         add_gene_tags_to_bam(
             assign_barcodes.out.chrom_bam_bc_bai
              // join on sample_id + chr
             .join(assign_genes.out.chrom_tsv_gene_assigns, by:[0, 1]))

         cluster_umis(add_gene_tags_to_bam.out.chrom_bam_bc_bai)

         // group by sample_id
         combine_chrom_bams(
             cluster_umis.out.bam_bc_bai.groupTuple())


         count_cell_gene_umi_reads(combine_chrom_bams.out.bam_fully_tagged)

         umi_gene_saturation(count_cell_gene_umi_reads.out.cell_umi_gene_tsv)

         construct_expression_matrix(combine_chrom_bams.out.bam_fully_tagged)

         process_expression_matrix(
            construct_expression_matrix.out.matrix_counts_tsv)

        preprocess_bams_for_stringtie(
            combine_chrom_bams.out.bam_fully_tagged
        )

        stringtie(
            gtf, 
            preprocess_bams_for_stringtie.out.flipped_bc_bams
            .flatMap{it ->
                l =[]
                // make [sample_id, bam] channel
                for (i=1; i<it[1].size(); i++) {
                    l.add( [it[0], it[1][i]] )
                }
                return l
            })

        build_transcript_matrix(
            stringtie.out.tmap
            // Collect all the tmap files by sample_id
            // Note: KeepHeader: false does not work here
            .collectFile()
            .map {it ->
                tuple(it.getBaseName(), it)}
            .join(process_expression_matrix.out.matrix_processed_tsv))

         umap_reduce_expression_matrix(
            process_expression_matrix.out.matrix_processed_tsv)

         umap_plot_total_umis(
             umap_reduce_expression_matrix.out.matrix_umap_tsv
             .join(process_expression_matrix.out.matrix_processed_tsv))

         genes_to_plot = Channel.fromPath(umap_genes)
             .splitCsv()
        
         umap_plot_genes(
             umap_reduce_expression_matrix.out.matrix_umap_tsv
             .join(process_expression_matrix.out.matrix_processed_tsv)
             .combine(genes_to_plot))

         umap_plot_mito_genes(
            umap_reduce_expression_matrix.out.matrix_umap_tsv
            .join(process_expression_matrix.out.matrix_mito_tsv))

     emit:
         umap_plots = umap_plot_genes.out
            .concat(
                umap_plot_total_umis.out.matrix_umap_plot_total
            )
            .groupTuple()
            .flatMap({it -> 
                l = []
                for (i=0; i<it[1].size(); i++)
                    if(l){
                        l.add(tuple(it[0], it[1][i]))
                    }
                return l
                })
                .concat(umap_plot_total_umis.out)
                .concat(umap_plot_mito_genes.out)
                .groupTuple()
                
        results = umi_gene_saturation.out
             .join(construct_expression_matrix.out)
             .join(process_expression_matrix.out.matrix_processed_tsv)
             .join(process_expression_matrix.out.matrix_mito_tsv)
             .join(generate_whitelist.out.whitelist)
             .join(generate_whitelist.out.kneeplot)
             .join(combine_chrom_bams.out.bam_fully_tagged)
             .join(extract_barcodes.out.barcode_counts)
             .join(umap_reduce_expression_matrix.out.matrix_umap_tsv)
}
