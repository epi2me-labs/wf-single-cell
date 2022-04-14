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


process summariseReads {
    // concatenate fastq and fastq.gz in a dir

    label "wftemplate"
    cpus 1
    input:
        tuple path(directory), val(sample_id), val(type)
    output:
        path "${sample_id}.stats"
    shell:
    """
    fastcat -s ${sample_id} -r ${sample_id}.stats -x ${directory} > /dev/null
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
    main:
        summary = summariseReads(reads)
        software_versions = getVersions()
        workflow_params = getParams()
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
    samples = fastq_ingress(
        params.fastq, params.out_dir, params.sample, params.sample_sheet, params.sanitize_fastq)

    pipeline(samples)
    output(pipeline.out.results)
    end_ping(pipeline.out.telemetry)
}
