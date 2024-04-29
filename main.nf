#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList;
nextflow.enable.dsl = 2

include { fastq_ingress } from './lib/ingress'
include { preprocess } from './subworkflows/preprocess'
include { process_bams } from './subworkflows/process_bams'

OPTIONAL_FILE = file("$projectDir/data/OPTIONAL_FILE")

process getVersions {
    label "singlecell"
    cpus 1
    memory "1 GB"
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import parasail; print(f'parasail,{parasail.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import rapidfuzz; print(f'rapidfuzz,{rapidfuzz.__version__}')" >> versions.txt
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
    memory "1 GB"
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
    label "wf_common"
    cpus 1
    memory "32 GB"
    input:
        val metadata
        path 'versions'
        path 'params.csv'
        path stats, stageAs: "stats_*"
        path 'survival.tsv'
        path 'wf_summary.tsv'
        path umap_dirs
        path images
        path umap_genes
        val wf_version

    output:
        path "wf-single-cell-*.html"
    script:
        String report_name = "wf-single-cell-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
    """
    echo '${metadata}' > metadata.json
    workflow-glue report \
        $report_name \
        --stats $stats \
        --params params.csv \
        --versions versions \
        --survival survival.tsv \
        --wf_summary wf_summary.tsv \
        --umap_dirs $umap_dirs \
        --images $images \
        --umap_genes $umap_genes \
        --metadata metadata.json \
        --wf_version $wf_version \
        --metadata metadata.json
    """
}


// See https://github.com/nextflow-io/nextflow/issues/1636
// This is the only way to publish files from a workflow whilst
// decoupling the publish from the process steps.
process output {
    label "singlecell"
    cpus 1
    memory "1 GB"
    // // publish inputs to output directory
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*.{bam,bai}",
        saveAs: { filename -> "${meta.alias}/bams/$filename" }
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*umap*.{tsv,png}",
        saveAs: { filename -> "${meta.alias}/umap/$filename" }
    publishDir "${params.out_dir}", mode: 'copy', 
        pattern: "*{images,counts,gene_expression,transcript_expression,kneeplot,saturation,config,tags,whitelist,transcriptome,annotation}*",
        saveAs: { filename -> "${meta.alias}/$filename" }

    input:
        tuple val(meta),
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
    cpus 1
    memory "1 GB"
    publishDir "${params.out_dir}", mode: 'copy', pattern: "*"
    input:
        path fname
    output:
        path fname
    """
    echo "Writing output files."
    """
}


process parse_kit_metadata {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path 'sample_ids'
        path sc_sample_sheet
        path kit_config, stageAs: 'kit_config.csv'
    output:
        path "merged.csv"
    script:
    if (sc_sample_sheet.name != "OPTIONAL_FILE"){
        """
        workflow-glue parse_kit_metadata from_sheet \
            --user_config ${sc_sample_sheet} \
            --kit_config kit_config.csv \
            --sample_ids sample_ids \
            --output merged.csv
        """
    }else{
        """
        workflow-glue parse_kit_metadata from_cli \
            --kit_config kit_config.csv \
            --kit_name "$params.kit_name" \
            --kit_version $params.kit_version \
            --expected_cells $params.expected_cells \
            --sample_ids $sample_ids \
            --output merged.csv
        """
    }

}


process prepare_report_data {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        tuple val(meta),
              path('read_tags'),
              path('config_stats'),
              path('white_list'),
              path('gene_mean_expression.tsv'),
              path('transcript_mean_expression.tsv'),
              path('mitochondrial_expression.tsv'),
              path(umaps)
    output:
        path 'survival_data.tsv',
            emit: survival
        path 'sample_summary.tsv',
            emit: summary
        path "${meta.alias}_umap",
            emit: umap_dir

    script:
        opt_umap = umaps.name != 'OPTIONAL_FILE'
        String hist_dir = "histogram_stats/${meta.alias}"
    """
    workflow-glue prepare_report_data \
        --read_tags read_tags \
        --config_stats config_stats \
        --white_list white_list \
        --sample_id ${meta.alias} \
        --summary_out sample_summary.tsv \
        --survival_out survival_data.tsv
    
    umd=${meta.alias}_umap
    mkdir \$umd

    if [ "$opt_umap" = true ]; then
        echo "Adding umap data to sample directory"
        # Add data required for umap plotting into sample directory
        mv *umap*.tsv \$umd
        mv gene_mean_expression.tsv \$umd
        mv transcript_mean_expression.tsv \$umd
        mv mitochondrial_expression.tsv \$umd
    else
        touch "\$umd"/OPTIONAL_FILE
    fi
    """
}


// workflow module
workflow pipeline {
    take:
        chunks
        ref_genome_dir
        umap_genes
    main:
        // throw an exception for deprecated conda users
        if (workflow.profile.contains("conda")) {
            throw new Exception(
                "Sorry, this workflow is not compatible with --profile conda," + 
                "please use --profile standard (Docker) " +
                "or --profile singularity.")
        }
        ref_genome_fasta = file("${params.ref_genome_dir}/fasta/genome.fa", checkIfExists: true)
        ref_genome_idx = file("${params.ref_genome_dir}/fasta/genome.fa.fai", checkIfExists: true)
        ref_genes_gtf = file("${params.ref_genome_dir}/genes/genes.gtf", checkIfExists: true)
        software_versions = getVersions()
        workflow_params = getParams()

        bc_longlist_dir = file("${projectDir}/data", checkIfExists: true)

        preprocess(
            chunks.map{meta, fastq, stats -> [meta, fastq]},
            bc_longlist_dir,
            ref_genome_fasta,
            ref_genome_idx,
            ref_genes_gtf)

        process_bams(
            preprocess.out.bam_sort,
            preprocess.out.read_tags,
            preprocess.out.high_qual_bc_counts.groupTuple(),
            ref_genes_gtf,
            ref_genome_fasta,
            ref_genome_idx)

        prepare_report_data(
            process_bams.out.final_read_tags
            .join(preprocess.out.config_stats)
            .join(process_bams.out.white_list)
            .join(process_bams.out.gene_mean_expression)
            .join(process_bams.out.transcript_mean_expression)
            .join(process_bams.out.mitochondrial_expression)
            .join(process_bams.out.umap_matrices))

        // Get the metadata and stats for the report
        chunks
            .groupTuple()
            .multiMap{ meta, chunk, stats ->
                meta: meta
                stats: stats[0]
            }.set { for_report }
        metadata = for_report.meta.collect()
        stats = for_report.stats.collect()

        makeReport(
            metadata,
            software_versions,
            workflow_params,
            stats,
            prepare_report_data.out.survival
                .collectFile(keepHeader:true),
            prepare_report_data.out.summary
                .collectFile(keepHeader:true),
            prepare_report_data.out.umap_dir.collect(),
            process_bams.out.plots,
            umap_genes,
            workflow.manifest.version)
    emit:
        results = process_bams.out.results
        config_stats = preprocess.out.config_stats
        report = makeReport.out
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)

    if (params.umap_plot_genes){
        umap_genes = file(params.umap_plot_genes, checkIfExists: true)
    }else{
        umap_genes = file("${projectDir}/umap_plot_genes.csv", checkIfExists: true)
    }

    if (params.kit_config){
        kit_configs_file = file(params.kit_config, checkIfExists: true)
    }else{
        kit_configs_file = file("${projectDir}/kit_configs.csv", checkIfExists: true)
    }

    fastq = file(params.fastq, type: "file")

    samples = fastq_ingress([
            "input":params.fastq,
            "sample":params.sample,
            "sample_sheet":params.sample_sheet,
            "fastq_chunk": params.fastq_chunk,
            "stats": true,
            "per_read_stats": false])


    if (!params.single_cell_sample_sheet) {

        sc_sample_sheet = file("$projectDir/data/OPTIONAL_FILE")
    } else {
        // Read single_cell_sample_sheet
        sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
    }

    fastqingress_ids = samples.map {meta, file, stats -> meta.alias }.unique().collectFile(newLine: true)
    // Get [sample_id, kit_meta]
    kit_meta = parse_kit_metadata(fastqingress_ids, sc_sample_sheet, kit_configs_file)
        .splitCsv(header:true)
        .map {it -> [it['sample_id'], it]}
    // Merge the kit metadata onto the sample metadata
    sample_and_kit_meta = kit_meta
        .cross(samples
            // Put sample_id as first element for join
            .map {meta, chunk, stats -> [meta.alias, meta, chunk, stats]})
        // Extract the joined sample and kit info from the cross results
        .map {kit, sample -> [ sample[1] + kit[1], sample[2], sample[3]]}
        // we never need the chunk index for merging items so discard it
        .map {meta, chunk, stats ->
            def new_meta = meta.clone()
            new_meta.remove('group_index')
            [new_meta, chunk, stats]}

    pipeline(
        sample_and_kit_meta,
        ref_genome_dir,
        umap_genes)

    output(pipeline.out.results.flatMap({it ->
        // Convert [meta, file, file, ..]
        // to      [[meta, file], [meta, file], ...]
        l = [];
            for (i=1; i<it.size(); i++) {
                l.add(tuple(it[0], it[i]))
            }
            return l
        }).concat(pipeline.out.config_stats)
    )

    output_report(pipeline.out.report)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}
