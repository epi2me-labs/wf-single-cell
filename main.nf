#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList;
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
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
    cache false
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
        path expression_dirs
        path images
        path umap_genes
        val wf_version
        path 'bam_stats.tsv'
        path visium_coords

    output:
        path "wf-single-cell-*.html"
    script:
        String report_name = "wf-single-cell-report.html"
        String metadata = new JsonBuilder(metadata).toPrettyString()
        def visium_opt = visium_coords.fileName.name != OPTIONAL_FILE.name ? '--visium_spatial_coords ' + visium_coords : ""
        String q_filtered = params.min_read_qual ? "--q_filtered": ""
    """
    echo '${metadata}' > metadata.json
    workflow-glue report \
        $report_name \
        --stats $stats \
        --params params.csv \
        --versions versions \
        --survival survival.tsv \
        --expr_dirs $expression_dirs \
        --images $images \
        --umap_genes $umap_genes \
        --metadata metadata.json \
        --wf_version $wf_version \
        --metadata metadata.json \
        --bam_stats bam_stats.tsv \
        $q_filtered \
        $visium_opt
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
        // A visium barcode is a tissue coordinate not a cell, so we don't need expected cells.
        if (params.kit.split(':')[0] != "visium" & params.expected_cells == null ){
            throw new Exception("expected_cells should be provided for 10x kits other than Visium")
        }
        """
        workflow-glue parse_kit_metadata from_cli \
            --kit_config kit_config.csv \
            --kit "$params.kit" \
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
              path('raw_gene_expression'),
              path('gene_mean_expression.tsv'),
              path('transcript_mean_expression.tsv'),
              path('mitochondrial_expression.tsv'),
              path(umaps),
              path('bamstats/bam_stats*.tsv')
        path("genes_of_interest.tsv")
    output:
        // sample_id column added to survival.tsv and bm_stats.tsv no need for meta
        path "survival.tsv", emit: survival
        path "bam_stats.tsv", emit: bam_stats
        path "${meta.alias}_expression", emit: expression_dir

    script:
        opt_umap = umaps.name != 'OPTIONAL_FILE'
        String hist_dir = "histogram_stats/${meta.alias}"
    """
    # Make a directory to stick some expression related files per sample
    expression_dir="${meta.alias}_expression"
    mkdir \$expression_dir
    echo \$expression_dir
    workflow-glue prepare_report_data \
        "${meta.alias}" adapter_stats bamstats expression_stats \
        white_list.txt survival.tsv bam_stats.tsv raw_gene_expression \
        genes_of_interest.tsv ${meta.n_seqs}

    if [ "$opt_umap" = true ]; then
        echo "Adding umap data to sample directory"
        # Add data required for umap plotting into sample directory
        mv *umap*.tsv \$expression_dir
        mv gene_mean_expression.tsv \$expression_dir
        mv transcript_mean_expression.tsv \$expression_dir
        mv mitochondrial_expression.tsv \$expression_dir
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
        genes_of_interest
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
        if (params.kit == 'visium:v1'){
            visium_coords = file("${bc_longlist_dir}/visium-v1_coordinates.txt", checkIfExists: true)
        }
        else {
             visium_coords = OPTIONAL_FILE
        }

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
            .join(process_bams.out.raw_gene_expression)
            .join(process_bams.out.gene_mean_expression)
            .join(process_bams.out.transcript_mean_expression)
            .join(process_bams.out.mitochondrial_expression)
            .join(process_bams.out.umap_matrices)
            .join(preprocess.out.bam_stats
                .groupTuple()),
            genes_of_interest
            )


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
            prepare_report_data.out.expression_dir.collect(),
            process_bams.out.plots,
            genes_of_interest,
            workflow.manifest.version,
            prepare_report_data.out.bam_stats
                .collectFile(keepHeader:true),
            visium_coords)
}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    Pinguscript.ping_start(nextflow, workflow, params)
    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)

    if (params.genes_of_interest){
        genes_of_interest = file(params.genes_of_interest, checkIfExists: true)
    }else{
        genes_of_interest = file("${projectDir}/data/genes_of_interest.csv", checkIfExists: true)
    }

    if (params.kit_config){
        kit_configs_file = file(params.kit_config, checkIfExists: true)
    }else{
        kit_configs_file = file("${projectDir}/kit_configs.csv", checkIfExists: true)
    }

    ArrayList fastcat_extra_args = []
    if (params.min_read_qual) {
        fastcat_extra_args << "-q $params.min_read_qual"
    }

    if (params.fastq) {
        samples = fastq_ingress([
                "input":params.fastq,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "fastq_chunk": params.fastq_chunk,
                "stats": true,
                "per_read_stats": false,
                "fastcat_extra_args": fastcat_extra_args.join(" ")])
    } else {
        samples = xam_ingress([
                "input":params.bam,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "fastq_chunk": params.fastq_chunk,
                "keep_unaligned": true,
                "return_fastq": true,
                "stats": true,
                "per_read_stats": false,
                "fastcat_extra_args": fastcat_extra_args.join(" ")])

    }

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
        genes_of_interest)
}

workflow.onComplete {
    Pinguscript.ping_complete(nextflow, workflow, params)
}
workflow.onError {
    Pinguscript.ping_error(nextflow, workflow, params)

}
