process get_ctat_data {
    label "wf_common"
    cpus 1
    memory "2 GB"
    storeDir {params.store_dir ? "${params.store_dir}/${name}" : null }
    input:
            val name
            val url
    output:
        path "${name}", emit: resource_dir
    script:
    """
    wget -qO- $url \
        | tar --no-same-owner -xzv --one-top-level=${name} --strip-component=1
    """
}

process find_fusions {
    /*
    Run ctat-LR-fusion to find fusion reads. 
    */
    label "ctat_lr_fusion"
    cpus params.wf.fusion_threads
    memory '16 GB'
    publishDir "${params.out_dir}/${meta.alias}/fusions", mode: 'copy', pattern: "${meta.alias}.ctat-LR-fusion.tar.gz"
    input:
        tuple val(meta),
              path("tagged.bam")
        path("ctat_reference_bundle")
    output:
        tuple val(meta),
              path("${meta.alias}.ctat-LR-fusion.tar.gz"),
              emit: gzipped_ctat_dir
        tuple val(meta),
              path(fusion_preds),
              emit: ctat_fusion_predictions
        stdout emit: stdout
    script:
        String ctat_outdir = "fusions"
        // Expected main output file
        // This will be present if no fusion candidates are verified (header only)
        // or absent if no fusion candidates are found
        fusion_preds = "${ctat_outdir}/ctat-LR-fusion.fusion_predictions.tsv"
        Integer threads = Math.max(2, task.cpus - 2)
    """
    ctat-LR-fusion \
        --LR_bam tagged.bam \
        --genome_lib_dir ./ctat_reference_bundle \
        --CPU ${threads} --vis --output ${ctat_outdir}

    if [ ! -f "${fusion_preds}" ]; then
        echo "No fusion candidates found for ${meta.alias}"
        # Create an empty file for expected output
        touch "${ctat_outdir}/ctat-LR-fusion.fusion_predictions.tsv"
    else
        n=\$(tail -n +2 "${fusion_preds}" | wc -l)
        if [ "\$n" -eq 0 ]; then
            echo "Fusion candidates found for ${meta.alias} but none passed filters"
        fi
    fi
    tar -czf "${meta.alias}.ctat-LR-fusion.tar.gz" ${ctat_outdir}
    """
}


process format_ctat_output {
    label "singlecell"
    cpus 1
    memory '2 GB'
    publishDir "${params.out_dir}/${meta.alias}/fusions", mode: 'copy', pattern: "${meta.alias}.ctat-LR-fusion.fusion_predictions_per*"
    input:
        tuple val(meta),
              path("ctat-LR-fusion.fusion_predictions.tsv"),
              path("read_summary_tags.tsv")
    output: 
        tuple val(meta),
              path("${meta.alias}.ctat-LR-fusion.fusion_predictions_per-read.tsv"),
              emit: read_summary
        tuple val(meta),
              path("${meta.alias}.ctat-LR-fusion.fusion_predictions_per-fusion.tsv"),
              emit: fusion_summary
        tuple val(meta),
              path("${meta.alias}.ctat-LR-fusion.fusion_summary.tsv"),
              emit: cell_summary
    script:
    """
    workflow-glue format_ctat_output \
        ctat-LR-fusion.fusion_predictions.tsv \
        read_summary_tags.tsv \
        "${meta.alias}.ctat-LR-fusion.fusion_predictions_per-read.tsv" \
        "${meta.alias}.ctat-LR-fusion.fusion_predictions_per-fusion.tsv" \
        "${meta.alias}.ctat-LR-fusion.fusion_summary.tsv" \
        ${meta.alias}
    """
}



workflow ctat_lr_fusion {
    take: 
        tagged_bam_and_summary
        ctat_reference_bundle
    main:
        find_fusions(
            tagged_bam_and_summary.map{meta, bam, _bai, _read_summary -> [meta, bam]}, 
            ctat_reference_bundle)
        
        find_fusions.out.stdout.map {stdout -> 
            if (stdout) {
                log.warn(stdout)
            }
        }
        
        format_ctat_output(
            find_fusions.out.ctat_fusion_predictions
                .join(tagged_bam_and_summary
                    .map {meta, _bam, _bai, read_summary -> [meta, read_summary]})
        )

    emit:
        read_summary = format_ctat_output.out.read_summary
        fusion_summary = format_ctat_output.out.fusion_summary
        cell_summary = format_ctat_output.out.cell_summary

}