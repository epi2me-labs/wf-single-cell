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
include { stranding } from './subworkflows/stranding'
include { stranding } from './subworkflows/align'


process summariseAndCatReads {
    // concatenate fastq and fastq.gz in a dir

    label "wftemplate"
    cpus 1
    input:
        tuple path(directory), val(meta)
    output:
        tuple val("${meta.sample_id}"), path("${meta.sample_id}.stats"), emit: stats
        tuple val("${meta.sample_id}"), path("${meta.sample_id}.fastq"), emit: fastq
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} > ${sample_id}.fastq
    """
}

// process get_kits {
//     // Add kit to the fastq tuple
//     cous 1
//     input:
//         tuple val(sample_id), path(fastq)
//         val sc_sample_sheet
//     output:
//          tuple val(sample_id), path(fastq), val(kit)
//     """
//     #!/usr/bin/env python
//     import pandas as pd

//     df = pd.read_csv($sc_sample_sheet, index_col=0)
//     kit_name = df.loc[$sample_id, 'kit_name']
//     kit_version = df.loc[$sample_id, 'kit_version']
//     """
// }

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
        sc_sample_sheet
        ref_genome_dir
    main:
        inputs = Channel.fromPath(sc_sample_sheet)
                    .splitCsv(header:true)
                    .map { row -> tuple(
                              row.run_id, 
                              row.kit_name, 
                              row.kit_version, 
                              file(row.path))}
        
        // Sockeye
        stranding(
            inputs,
            sc_sample_sheet)
        
        // 10x reference downloads have known names
        REF_GENOME_FASTA = params.REF_GENOME_DIR / "fasta/genome.fa"
        REF_GENES_GTF = params.REF_GENOME_DIR / "genes/genes.gtf"
        
        align(
            stranding.out.STRANDED_FQ,
            REF_GENOME_FASTA,
            REF_GENES_GTF
        )
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {
    sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)

    pipeline(sc_sample_sheet, ref_genome_dir)
}
