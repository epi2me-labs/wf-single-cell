// Merge TSVs and sum the specified column
// Currently support only headerless inputs and summing of the second column
process merge_and_publish_tsv {
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    label "wf_common"
    cpus 1
    memory "2GB"
    input:
        tuple val(meta),
              path("inputs/input*.tsv")
        val(output_fname)
    output:
        tuple val(meta),
              path("${meta.alias}.${output_fname}")
    script:
    """
    find inputs -name "*.tsv" \
        -exec cat {} + \
        | csvtk -t summary -H -f 2:sum -g 1 \
        > "${meta.alias}.${output_fname}"
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

process call_paftools {
    label "singlecell"
    memory "2 GB"
    cpus 1
    input:
        path "ref_genes.gtf"
    output:
        path "ref_genes.bed", emit: ref_genes_bed
    script:
    """
    paftools.js gff2bed -j ref_genes.gtf > ref_genes.bed
    """
}

process cat_tags_by_chrom {
    // Merge per-chunk tags to create per-chromosome tags
    label "wf_common"
    cpus params.threads
    memory "8 GB"
    input:
        tuple val(meta),
              path('tags/*tags.tsv')
    output:
         tuple val(meta),
              path("chr_tags/*"),
              emit: merged_tags

    script:
    """
    mkdir chr_tags
    # Find the chr column number
    files=(tags/*)
    chr_col=\$(awk -v RS='\t' '/chr/{print NR; exit}' "\${files[0]}")

    # merge the tags TSVs, keep header from first file and split entries by chromosome
    awk -F'\t' -v chr_col=\$chr_col 'FNR==1{hdr=\$0; next} \
    {if (!seen[\$chr_col]++) \
        print hdr>"chr_tags/"\$chr_col".tsv"; \
        print>"chr_tags/"\$chr_col".tsv"}' tags/*

    # Sort by corrected barcode to allow chunked reading later.
    for file in chr_tags/*.tsv;
    do
        csvtk sort -t --keys CB "\$file" -o tmp.tsv;
        mv tmp.tsv "\$file";
    done
    """
}