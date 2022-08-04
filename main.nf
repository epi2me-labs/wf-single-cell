#!/usr/bin/env nextflow

// Developer notes
//
// This template workflow provides a basic structure to copy in order
// to create a new workflow. Current recommended pratices are:
//     i) create a simple command-line interface.
//    ii) include an abstract workflow scope named "pipeline" to be used
//        in a module fashion
//   iii) a second concreate, but anonymous, workflow scope to be used
//        as an entry point when using this workflow in isolation.

import groovy.json.JsonBuilder
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { start_ping; end_ping } from './lib/ping'
include { stranding } from './subworkflows/stranding'


process summariseAndCatReads {
    // concatenate fastq and fastq.gz in a dir

    label "wftemplate"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val(meta.sample_id}, path(${meta.sample_id}.stats), emit: stats
        tuple val(meta.sample_id}, path(${meta.sample_id}.fastq), emit: fastq
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} > ${sample_id}.fastq
    """
}


process getVersions {
    label "wftemplate"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    """
}


process getParams {
    label "wftemplate"
    cpus 1
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}

process stranding {
    input:
        tuple val(sample_id), path(reads)
    output:

}



process makeReport {
    label "wftemplate"
    input:
        path "seqs.txt"
        path "versions/*"
        path "params.json"
    output:
        path "wf-template-*.html"
    script:
        report_name = "wf-template-" + params.report_name + '.html'
    """
    report.py $report_name \
        --versions versions \
        seqs.txt \
        --params params.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    // publish inputs to output directory
    label "wftemplate"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        sc_sample_sheet
        ref_genome_dir
    main:
        summary = summariseReads(reads)
        software_versions = getVersions()
        workflow_params = getParams()
        
        // Sockeye
        stranding(reads)



        report = makeReport(summary, software_versions.collect(), workflow_params)
    emit:
        results = summary.concat(report)
        // TODO: use something more useful as telemetry
        telemetry = workflow_params
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    start_ping()

    sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)


    samples = fastq_ingress([
        "input":params.fastq,
        "sample":params.sample,
        "sample_sheet":params.sample_sheet,
        "sanitize": params.sanitize_fastq,
        "output":params.out_dir])

    pipeline(samples, sc_sample_sheet, ref_genome_dir)
    output(pipeline.out.results)
    end_ping(pipeline.out.telemetry)
}
