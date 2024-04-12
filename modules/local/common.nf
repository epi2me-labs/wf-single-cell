process construct_expression_matrix {
    label "singlecell"
    cpus  1
    memory = { 16.GB * task.attempt }
    maxRetries = 3
    errorStrategy = { task.exitStatus in 137..140 ? 'retry' : 'terminate' }
    input:
        tuple val(meta),
              path("read_tags.tsv")
        val(feature_type)
    output:
        tuple val(meta),
              val(feature_type),
              path("${meta.alias}.${feature_type}_expression.count.tsv"),
              emit: matrix_counts_tsv
        tuple val(meta),
              val(feature_type),
              path("${meta.alias}.${feature_type}_expression.processed.tsv"),
              emit: matrix_processed_tsv
        tuple val(meta),
              val(feature_type),
              path("${meta.alias}.${feature_type}_mean_per_cell_expression.tsv"),
              emit: mean_per_cell_expression
        tuple val(meta),
              path("${meta.alias}.${feature_type}_expression.mito.tsv"),
              emit: mito_expression_tsv,
              optional: true  // Only output for 'gene' feature_type
    """
    # Split the comma-separated mito prefixes
    IFS=","
    read -ra mito_prefixes <<< "${params.mito_prefix}"

    workflow-glue expression_matrix \
        read_tags.tsv \
        "${feature_type}" \
        "${meta.alias}.${feature_type}" \
        --min_features $params.matrix_min_genes \
        --min_cells $params.matrix_min_cells \
        --max_mito $params.matrix_max_mito \
        --mito_prefixes "\${mito_prefixes[@]}" \
        --norm_count $params.matrix_norm_count
    """
}
