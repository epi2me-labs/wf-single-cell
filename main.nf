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
include { align } from './subworkflows/align'
include { process_bams } from './subworkflows/process_bams'


process summariseCatChunkReads {
    // concatenate fastq and fastq.gz in a dir. 
    // Split into p parts where p is num threads

    label "singlecell"
    cpus 1
    input:
        tuple path(directory), val(meta)
        val check  // This will not exist if the sample_id check fails and will halt the pipleine.
    output:
        tuple val("${meta.sample_id}"), path("${meta.sample_id}.stats"), emit: stats
        tuple val("${meta.sample_id}"), path("chunks/*"), emit: fastq_chunks
    shell:
    """
    fastcat -s ${meta.sample_id} -r ${meta.sample_id}.stats -x ${directory} | \
        seqkit split2 -p ${params.max_threads} -O chunks -o ${meta.sample_id} -e .fastq
    """
}

process makeReport {
    label "singlecell"
    output:
        path "wf-single-cell-*.html"
    script:
        report_name = "wf-single-cell-" + params.report_name + '.html'
    """
    echo "<html><body>This report is a placeholder<br> \
        a detailed report will be available in the next release</body></html>" \
        > wf-single-cell-report.html
    """
}



// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    label "singlecell"
    // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*", 
        saveAs: { filename -> "${sample_id}/$filename" }
    label "isoforms"
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


process pack_images {
    label "singlecell"
    input:
        tuple val(sample_id),
              path(imgs)  
        output:
            tuple val(sample_id),
                path('umap_plots')
    """
    mkdir umap_plots
    for img in $imgs; do
        cp \$img umap_plots
    done;
    """
}

process check_sampleids{
    // Check that sample_ids gicven in the single_cell_sample_sheet are 
    // identical to the sample_ids of the fastq inputs
    input:
        path fastqingress_ids
        path sc_sample_sheet_ids
    output:
        // env check_sampleids_PASSED, emit: passed
        path 'diff', optional: true, emit: diff
    """
    #!/usr/bin/env python
    import pandas as pd
    import sys
    df_s = pd.read_csv("$fastqingress_ids", index_col=None)
    df_f = pd.read_csv("$sc_sample_sheet_ids", index_col=None)

    if set(df_s.iloc[:, 0].values) == set(df_f.iloc[:, 0].values):
        print('Success. The sample_ids are the same')
        open('diff', 'w').close()
    else:
        print("The smaples are different")
        sys.stdout.write('ksfhdskhjfsdkjhksjdaskjd')
    """
}


// workflow module
workflow pipeline {
    take:
        reads
        sc_sample_sheet
        ref_genome_dir
        umap_genes
        sample_kits
        sample_ids_check
    
    main:
        error_msg = ""
        ref_genome_fasta = file("${ref_genome_dir}/fasta/genome.fa", checkIfExists: true)
        ref_genes_gtf = file("${ref_genome_dir}/genes/genes.gtf", checkIfExists: true)
        ref_genome_idx = file("${ref_genome_fasta}.fai", checkIfExists: true)
        
        if (params.kit_config){
            kit_configs = file(params.kit_config, checkIfExists: true)
        }else{
            kit_configs = file("${projectDir}/kit_configs.csv", checkIfExists: true)
        }
        
        bc_longlist_dir = file("${projectDir}/data", checkIfExists: true)

        summariseCatChunkReads(reads, sample_ids_check)

        stranding(
            summariseCatChunkReads.out.fastq_chunks,
            sample_kits)

        align(
            stranding.out.stranded_fq,
            ref_genome_fasta,
            ref_genes_gtf,
            ref_genome_idx
        )

        process_bams(
            align.out.bam_sort,
            align.out.bam_sort_bai,
            sc_sample_sheet,
            kit_configs,
            ref_genes_gtf,
            umap_genes,
            bc_longlist_dir,
            sample_kits,
            ref_genome_fasta
        )
    emit:
        results = process_bams.out.results
        umap_plots = process_bams.out.umap_plots
        config_stats = stranding.out.config_stats
        
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    if (params.disable_ping == false) {
        try { 
            Pinguscript.ping_post(workflow, "start", "none", params.out_dir, params)
        } catch(RuntimeException e1) {
        }
    }

    sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)
    umap_genes = file(params.umap_plot_genes, checkIfExists: true)

    fastq = file(params.fastq, type: "file")

    reads = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "sanitize": params.sanitize_fastq,
            "output":params.out_dir])
    
    sample_kits = Channel.fromPath(sc_sample_sheet)
                    .splitCsv(header:true)
                    .map { row -> tuple(
                              row.sample_id, 
                              row.kit_name, 
                              row.kit_version,
                              row.exp_cells)}

    fastqingress_ids = reads.map{it -> it[1]['sample_id']}
    .collectFile(name: 'fastingress_read_ids.csv', newLine: true)
        
    sample_kit_ids = sample_kits.map{it -> it[0]}
        .collectFile(name: 'sc_sample_sheet_ids.csv', newLine: true)    

    check_sampleids(fastqingress_ids, sample_kit_ids)

    check_sampleids.out.ifEmpty{
        exit 1, 
        """
        The sample_ids in the single_cell_sample_sheet do not match those
        of the fastq inputs. Please see the README for instructions """.stripIndent()}

    // first() converts the queue channel to a value channel.
    pipeline(reads, sc_sample_sheet, ref_genome_dir, umap_genes, sample_kits,
        check_sampleids.out.first())

    pack_images(pipeline.out.umap_plots)
    
    output(pipeline.out.results.flatMap({it ->
        // Convert [sample_id, file, file, ..] 
        // to      [[sample_id, file], [sample_id, file], ...]
        l = [];
            for (i=1; i<it.size(); i++) {
                l.add(tuple(it[0], it[i]))
            }
            return l
        }).concat(pack_images.out,
            pipeline.out.config_stats)
    )
    // This is temporay until a detailed report is made
    makeReport()
    output_report(makeReport.out)

    if (params.disable_ping == false) {
    workflow.onComplete {
        try{
            Pinguscript.ping_post(workflow, "end", "none", params.out_dir, params)
        }catch(RuntimeException e1) {
        }
    }
    
        workflow.onError {
            try{
                Pinguscript.ping_post(workflow, "error", "$workflow.errorMessage", params.out_dir, params)
            }catch(RuntimeException e1) {
            }
        }

    }
    
    
}
