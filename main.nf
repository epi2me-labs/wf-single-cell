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
import nextflow.util.BlankSeparatedList;
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/fastqingress'
include { stranding } from './subworkflows/stranding'
include { align } from './subworkflows/align'
include { process_bams } from './subworkflows/process_bams'


process summariseCatChunkReads {
    // concatenate fastq and fastq.gz in a dir. 
    // Split into p parts where p is num threads

    label "singlecell"
    cpus 2
    input:
        tuple path(directory), val(meta)
    output:
        tuple val("${meta.sample_id}"), 
              path("${meta.sample_id}.stats"), 
              emit: stats
        tuple val("${meta.sample_id}"), 
              path("chunks/*"),
              emit: fastq_chunks
    
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} | \
        seqkit split2 -p ${params.max_threads} -O chunks -o ${meta.sample_id} -e .gz
    """
}

process getVersions {
    label "singlecell"
    cpus 1
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import sklearn; print(f'scikit-learn,{sklearn.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    stringtie --version | sed 's/^/stringtie,/' >> versions.txt
    gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "singlecell"
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
    label "singlecell"
    input:
        path 'versions'
        path 'params.csv'
        path 'read_stats.csv'
        path 'survival.tsv'
        path 'wf_summary.tsv'
        path images

    output:
        path "wf-single-cell-*.html"
    script:
        report_name = "wf-single-cell-" + params.report_name + '.html'
    """
    report.py \
        --read_stats read_stats.csv \
        --params params.csv \
        --versions versions \
        --survival survival.tsv \
        --wf_summary wf_summary.tsv \
        --output ${report_name} \
        --images ${images} \
    """
}



// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    label "singlecell"
    // // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.{bam,bai}",
        saveAs: { filename -> "${sample_id}/bams/$filename" }
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*umap*.{tsv,png}",
        saveAs: { filename -> "${sample_id}/umap/$filename" }
    publishDir "${params.out_dir}", mode: 'copy', 
        pattern: "*{counts,processed,kneeplot,saturation,config,tags}*",
        saveAs: { filename -> "${sample_id}/$filename" }

    input:
        tuple val(sample_id),
              path(fname)
    output:
        path fname
    """
    echo "Writing output files"
    """
}


process output_report {
    // publish inputs to output directory
    label "singlecell"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


process prepare_report_data {
    label "singlecell"
    input:
        tuple val(sample_id),
              path('read_tags'),
              path('config_stats'),
              path('white_list')
    output:
        path 'survival_data.tsv',
            emit: survival
        path 'sample_summary.tsv',
            emit: summary
    """
    prepare_report_data.py \
        --read_tags read_tags \
        --config_stats config_stats \
        --white_list white_list \
        --sample_id ${sample_id}

    """
}


// workflow module
workflow pipeline {
    take:
        reads
        ref_genome_dir
        umap_genes
        meta
    main:
        // throw an exception for deprecated conda users
        if (workflow.profile.contains("conda")) {
            throw new Exception(
                "Sorry, this workflow is not compatible with --profile conda," + 
                "please use --profile standard (Docker) " +
                "or --profile singularity.")
        }
        ref_genome_fasta = file("${ref_genome_dir}/fasta/genome.fa", checkIfExists: true)
        ref_genome_idx = file("${ref_genome_fasta}.fai", checkIfExists: true)
        ref_genes_gtf = file("${ref_genome_dir}/genes/genes.gtf", checkIfExists: true)
        software_versions = getVersions()
        workflow_params = getParams()
        
        bc_longlist_dir = file("${projectDir}/data", checkIfExists: true)
        
        summariseCatChunkReads(reads)

        stranding(
            summariseCatChunkReads.out.fastq_chunks, meta)

        align(
            stranding.out.stranded_fq,
            ref_genome_fasta,
            ref_genome_idx,
            ref_genes_gtf)

        process_bams(
            align.out.bam_sort,
            align.out.chr_beds,
            meta,
            ref_genes_gtf,
            umap_genes,
            bc_longlist_dir,
            ref_genome_fasta,
            ref_genome_idx)
        
        prepare_report_data(
            process_bams.out.read_tags
            .join(stranding.out.config_stats)
            .join(process_bams.out.white_list))
        
        makeReport(
            software_versions,
            workflow_params,
            summariseCatChunkReads.out.stats
                .map {it -> it[1]}
                .collectFile(keepHeader:true, name: 'stats'),
            prepare_report_data.out.survival
                .collectFile(keepHeader:true),
            prepare_report_data.out.summary
                .collectFile(keepHeader:true),
            process_bams.out.plots)
    emit:
        results = process_bams.out.results
        config_stats = stranding.out.config_stats
        report = makeReport.out
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
    }
    
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)
    umap_genes = file(params.umap_plot_genes, checkIfExists: true)

    if (params.kit_config){
        kit_configs_file = file(params.kit_config, checkIfExists: true)
    }else{
        kit_configs_file = file("${projectDir}/kit_configs.csv", checkIfExists: true)
    }

    fastq = file(params.fastq, type: "file")

    reads = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet])

    fastqingress_ids = reads.map{it -> it[1]['sample_id']}

    if (!params.single_cell_sample_sheet) {
        //build_a single_cell_sample_sheet channel applying the same kit values to each sample
        kit_string = Channel.value("${params.kit_name},${params.kit_version},${params.expected_cells}")
        .splitCsv() 
        sample_kits = fastqingress_ids
            .combine(kit_string)
            .map{ it -> 
                [it[1], it[2], 
                ["sample_id": it[0], 
                'kit_name': it[1], 
                'kit_version:': it[2], 
                'exp_cells': it[3]]] }
    } else {
        // Read single_cell_sample_sheet
        sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
        sample_kits = Channel.fromPath(sc_sample_sheet)
            .splitCsv(header:true)
            .map {it -> [it['kit_name'], it['kit_version'], it]}
    }

    kit_configs = Channel.fromPath(kit_configs_file)
        .splitCsv(header:true)
        .map {it -> [it['kit_name'], it['kit_version'], it]}

    // Do a check for mismatching sample_ids in the data and single cell sample sheet
    sample_kits.map {it -> it[2]['sample_id']}
    .join(fastqingress_ids, failOnMismatch:true)

    // Merge the kit info and user-supplied meta data on kit name and version
    sample_info = kit_configs.join(sample_kits, by: [0, 1])
        .map {it ->
            meta = it[2] + it[3] // Join the 2 meta maps  
            kit_name = meta['kit_name']
            kit_version = meta['kit_version']
            // and Get the appropriate cell barcode longlist based on the kit_name specified for this sample_id.
            switch(kit_name){
                case '3prime':
                    switch(kit_version){
                        case 'v2':
                            long_list = "737K-august-2016.txt.gz"
                            break
                        case 'v3':
                            long_list = "3M-february-2018.txt.gz"
                            break
                        default:
                            throw new Exception(
                                "Encountered an unexpected kit version for 3prime kit (v2 or v3): ${kit_version}")
                    }
                    break;
                case '5prime':
                    long_list = "737K-august-2016.txt.gz"
                    break
                case 'multiome':
                    long_list = "737K-arc-v1.txt.gz"
                    break
                default:
                    throw new Exception("Encountered an unexpected kit_name in samples.csv")
            }
            meta['bc_long_list'] = long_list

            [it[3]['sample_id'], meta]}

    pipeline(reads, ref_genome_dir, umap_genes, sample_info)



    output(pipeline.out.results.flatMap({it ->
        // Convert [sample_id, file, file, ..] 
        // to      [[sample_id, file], [sample_id, file], ...]
        l = [];
            for (i=1; i<it.size(); i++) {
                l.add(tuple(it[0], it[i]))
            }
            return l
        }).concat(pipeline.out.config_stats)
    )

    output_report(pipeline.out.report)
}

if (params.disable_ping == false) {
    workflow.onComplete {
        Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
    }

    workflow.onError {
        Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
    }

}
