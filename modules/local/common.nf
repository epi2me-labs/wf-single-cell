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