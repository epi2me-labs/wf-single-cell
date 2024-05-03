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
    publishDir "${params.out_dir}", mode: 'copy', pattern: "wf-single-cell-report.html"
    input:
        val metadata
        path 'versions'
        path 'params.csv'
        path stats, stageAs: "stats_*"
        path 'survival.tsv'
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
        --umap_dirs $umap_dirs \
        --images $images \
        --umap_genes $umap_genes \
        --metadata metadata.json \
        --wf_version $wf_version \
        --metadata metadata.json
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
              path('adapter_stats/stats*.json'),
              path('expression_stats/stats*.json'),
              path('white_list.txt'),
              path('gene_mean_expression.tsv'),
              path('transcript_mean_expression.tsv'),
              path('mitochondrial_expression.tsv'),
              path(umaps)
    output:
        path "survival.tsv", emit: survival  // not meta.alias here, breaks collectFile()
        path "${meta.alias}_umap", emit: umap_dir

    script:
        opt_umap = umaps.name != 'OPTIONAL_FILE'
        String hist_dir = "histogram_stats/${meta.alias}"
    """
    workflow-glue prepare_report_data \
        "${meta.alias}" adapter_stats expression_stats white_list.txt survival.tsv
    
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
            preprocess.out.adapter_summary.groupTuple()
            .join(process_bams.out.expression_stats
                .groupTuple()
                .map{meta, chrs, stats -> [meta, stats]})
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

        // note the cheeky little .collectFile() here to concatenate the
        // read survival stats from different samples into a single file
        makeReport(
            metadata,
            software_versions,
            workflow_params,
            stats,
            prepare_report_data.out.survival
                .collectFile(keepHeader:true),
            prepare_report_data.out.umap_dir.collect(),
            process_bams.out.plots,
            umap_genes,
            workflow.manifest.version)
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
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}
