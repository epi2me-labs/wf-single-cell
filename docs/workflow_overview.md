# wf-single-cell workflow overview

The purpose of this document is to provide a detailed 
overview of the `wf-single-cell` workflow. 
By the end of this document, you should have a good understanding of the main workflow steps, 
and be equipped with the knowledge you need to use it effectively.

Reads generated from 10x Genomics libraries contain cell barcodes and unique molecular identifiers (UMIs),
which enable the assignment of the read to a cell and a unique RNA molecule. 
The main purpose of this workflow is to extract and correct barcodes and UMIs derived from 10x Genomics libraries 
that have been sequenced with Oxford Nanopore Technologies (ONT) instruments. Gene and transcript assignment 
is then performed on the corrected reads, resulting in gene/transcript x cell count matrices, which can be used in 
downstream applications. 

Workflow parameters specified in this document are shown in the command line format:
`--param_name value`. Parameters can also be defined in a Nextflow config file like so: `param_name = value`. The former is used in this document. 
For more information see the [Configuration section of the Nextflow documentation](https://www.nextflow.io/docs/latest/config.html).

 
## The main stages of the workflow
* Stranding and identification of full-length reads
* Aligning reads to genome
* Extract cell barcodes and UMIs
* Barcode correction
* Gene and transcript assignment
* UMI correction
* Make expression matrices
* Tagging BAM files
* Calculate library saturation
* Make UMAP plots

## Stranding and identification of full-length reads
Reads derived from a 10x Genonomics library will ideally be flanked by two different adapter
sequences. These reads are more likely to represent full length mRNA molecules,
although that isn't guaranteed as some cDNAs may have been created by internal
priming or represent other cDNA synthesis artifacts. See this [10x Genomics Technical Note](https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf) for more information.


This stage of the workflow assigns an adapter configuration to each read, 
as well as determining the strand of the read and orientating the reads into the same orientation. 
The following schematic shows an example read structure from 10x Genomics 3&#8242; library .

<figure>
<img src="images/3prime_read.png", alt="10x read structure">
<figcaption>Fig.1 Read structure for 10x 3prime kit reads</figcaption>
</figure>

Adapter are located within the reads using [vsearch](https://github.com/torognes/vsearch) (`Read1` and `TSO` in Fig.1 in the case of the 3prime kit).
The table below details the adapter sequences for each the 10x Genomics kits, along with links to the relevant user guides.

<br> 

| Kit      | adapter1               | adapter2 | 10x user guide |
|----------|------------------------|----------|----------------|
| 3&#8242;   | CTACACGACGCTCTTCCGATCT | ATGTACTCTGCGTTGATACCACTGCTT | [3' kit](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf) |
| multiome | CTACACGACGCTCTTCCGATCT | ATGTACTCTGCGTTGATACCACTGCTT | [5'  kit](https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000331_ChromiumNextGEMSingleCell5-v2_UserGuide_RevE.pdf) |
| 5&#8242;   | CTACACGACGCTCTTCCGATCT | GTACTCTGCGTTGATACCACTGCTT | [multiome kit](https://teichlab.github.io/scg_lib_structs/data/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevB.pdf) |

Concatenated reads are identified by adapter configuration and are split into individual subreads and 
reorientated if required.
The following table details the various configurations and the actions taken for each.


| configuration   | action                                                    |
|-----------------|-----------------------------------------------------------|
| Full length     | Reads are trimmed from each side and oriented             |
| Single adapter1 | Reads are oriented and trimmed from adapter1 end only     |
| Single adapter2 | Reads are oriented and trimmed from the adapter2 end only |
| double adapter1 | Reads are trimmed from both sides                         |
| double adapter2 | Reads are trimmed from both sides                         |
| other           | No valid adapters found; not used in further analysis     | 


## Aligning reads to genome
The next stage is to align the preprocessed reads to the reference genome. This enables gene and transcript
read assignment in downstream steps.

The stranded fastq reads are mapped to the reference genome using minimap2.  
`--resources_mm2_max_threads int` controls the threads given to an alignment process.
Other optional parameters can be supplied to minimap2 using `--resources_mm2_flags` (for example `--resources_mm2_flags "-I 16GB"`).

## Extract cell barcodes and UMIs
The next step is to extract 10x Genomics barcodes and UMI sequences from the stranded and trimmed reads.

In order to do this, the first 100bp of each read are aligned to a reference probe using [parasail](https://github.com/jeffdaily/parasail). This probe contains a suffix of the adapter1 sequence, some ambiguities ("Ns") representing the barcode and UMI, and a polyT tract.

<figure>
<img src="images/probe.png" alt="probe image">
<figcaption>Fig.2 Schematic of a 3&#8242; v3 probe aligned to the read prefix</figcaption>
</figure>

In this case, the first 16bp of the read that aligns to the N probe region is where the barcode is expected to be, so
this sequence is assigned as `uncorrected barcode`.
The following 12 bp are where the UMI is expected, and is assigned to `uncorrected umi`.

The probe sequence will vary depending on the size of the expected barcode and UMI, which can be found in the 
following table.

| kit      |version|barcode length| UMI length |
|----------|---|---|------------|
| 3&#8242;   |v2|16| 10         |
| 3&#8242;   |v3|16| 12         |
| 5&#8242;   |v1|16| 10         |
| multiome |v1|16| 12         |

Workflow options:
- The size of the adapter1 suffix can be specified with: `--barcode_adapter1_suff_length`

## Barcode correction
The aim of this stage is to correct errors present in the previously extracted barcodes.
10x Genomics barcodes are not random; all possible barcode sequences can be found in a whitelist of known barcodes. 
The approprate whitelist for each kit and version is chosen automatically by the workflow.
The whitelist is used to generate a shortlist of high quality barcodes present in the sample, which is then used to 
correct barcode errors.

The correction proceeds as follows:
* A shortlist of all high quality barcodes present in our sample is generated. 
This is done by adding an uncorrected barcode to the shortlist if it:
  * has a min quality > `--barcode_min_quality` (default 15)

In each cell library there are expected to be some low quality cells and empty droplets that can be identified by their low number of reads.
To remove these cells, the shortlist is filtered with a quantile based threshold.
This threshold is determined by ranking the cells by read count and taking the top n cells (n = `--expected_cells`).
The read count 95th percentile / 20 is the threshold used. This threshold can be visualised in the knee plots generated by the workflow.

For uncorrected barcodes not present in the shortlist, they are cross-referenced against the shortlist, and are assigned 
a barocode from this list if the candidate barcode meets the following criteria:
* the query and closest-matched shortlist barcode have an edit distance <= 2
* The edit distance difference between the top and hit and second top hit in the shortlist >= 2

workflow options:
* `--barcode_max_ed`: Max edit distance between the uncorrected barcode and the matching whitelist barcode (default 2)".
* `--barcode_min_ed_diff`: Min difference in edit distance between (1) the uncorrected barcode vs top hit and (2) uncorrected barcode vs runner-up hit (default 2)
* `--expected_cells`: Number of expected cells. Enter the number of expected cells in the sample


## Gene and transcript assignment
Now that barcodes and UMIs have been assigned, the next stage is to assign a gene and transcript to each read. 

The first step is to create a transcriptome using stringtie2 (using long read mode `-L`).
The minimum required transcripts for a transcript to be called is 2. This can be changed, along with other stringtie settings using the 
`--stringtie_opts` options (default "-c 2"). The resulting query transcriptome is then annotated with the supplied reference 
genome annotation using [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), annotating query transcripts with reference transcript and gene metadata.

The original reads are then mapped to this transcriptome, producing one or more alignments per read.

### Transcript assignment
To assign a read to a transcript, we follow a similar procedure to that used by [FLAMES](https://github.com/LuyiTian/FLAMES).
* Transcripts that map to intronic query transcripts (those mapping to query transcripts with gffcompare classes ['i', 'y', 'p', 's'])
are set to unknown (-). This does not affect the gene assignment.
* Any read that uniquely maps to the annotated transcriptome is assigned that transcript.
* If the read does not uniquely map to a transcript, use the following procedure:
  * Order alignments by alignment score (AS)
  * If the top alignment has a better AS or query coverage than the second, and
  has a minimum transcript coverage of 0.4, assign to this transcript, else set as unknown (-)
  * If the AS is the same between the 1st and 2nd alignment 
    * Choose the alignment that has the highest transcript coverage, if the transcript coverage is > 0.8 

### Gene assignment
For each read, the alignment with the highest genomic MAPQ score is used for gene assignemnt if this score > `--gene_assigns_minqv` (default 30).
The gene ID is then derived from the transcript assigned in the transcript assignment step.

## UMI correction 
The next step is to correct UMI errors. The previous barcode correction step leveraged a
whitelist of all potential barcodes, drastically narrowing the search space. However, UMIs are random sequences and so there is no way of knowing beforehand
which UMIs will be in our library. To reduce the UMI search space, the reads are initially clustered by assigned gene, or if
no gene is assigned, by genomic interval. This also has the effect of reducing the likelihood of UMI collisions which occur due to the number
of total reads in a library usually being significantly higher than the total possible number of UMIS (16.7 million combinations for a 12bp UMI).

UMI clusters are then generated from these gene-clustered reads using [UMI-tools](https://github.com/CGATOxford/UMI-tools) with 
a slightly modified version of the [directional method](https://umi-tools.readthedocs.io/en/latest/the_methods.html).
For this clustering, Levenshtein distance threshold is used (default 2) instead of hamming distance, in order to account for UMI indels. 
The selected node of each cluster is the one with the highest number of reads and is used to assign a `corrected_umi` to the rest of the reads in that cluster.

## Make expression matrices
The gene x cell and transcript x cell expression matrices are the main outputs of the workflow and can be used in further analysis of 
the single cell experiment. 

The expression matrices are generated by collapsing the corrected UMIs (counting each unique UMI once) and summing
the number per feature/cell to give a feature/cell expression matrix (`*expression.counts.tsv`). 

The expression count matrices are further processed in the following way to give gene x cell processed matrices 
(`*expression.processed.tsv`):
* Cells are dropped that contain less than `--matrix_min_genes` genes or transcripts (default 200)
* Genes are dropped which are present in fewer than `--matrix_min_cells` (default 3)
* Cells where mitochondrial genes make up more than `--matrix_max_mito` (default 20%) are dropped
* Counts are normalized to `--matrix_norm_count` (default 10,000) reads/cell
* Normalized counts are finally log10 transformed


## Tagging bam files
BAM files generated from aligning reads to the reference genome are now tagged with the 
following information from the workflow. By default, the BAMs are output per chromosome, but can be concatenated
into a single file per sampe using `--merge_bam`
  - CB: corrected cell barcode sequence
  - CR: uncorrected cell barcode sequence
  - CY: Phred quality scores of the uncorrected cell barcode sequence
  - UB: corrected UMI sequence
  - UR: uncorrected UMI sequence
  - UY: Phred quality scores of the uncorrected UMI sequence


## Calculate library saturation
Sequencing saturation is an estimate of the library complexity that has been captured in a sequencing run.
As read depth per cell increases, the number of genes and distinct UMIs identified will increase at a rate that is 
dependent on the complexity of the input library. 

Reads which have been assigned a corrected barcode, corrected UMI and gene
are subsampled by varying degrees and the resulting median reads/cell is plotted against either:
- median genes per cell: this gives an indication of the gene complexity captured and whether increasing read depth would lead 
 to significantly more genes being identified.
- median UMIs per cell: this gives an indication of the UMI complexity captured. If this plot plateaus out, it indicates that 
    PCR duplicates represent a high proportion of the reads at higher sequencing depth.
- Sequencing saturation:  calculated as `1 - (number of unique UMIs across all cells / number of reads across all cells)`. Values near 0 indicate that
 very few duplicate UMIs are being identified, whereas values nearer to 1 indicate a higher proportion of duplicate UMIs are being captured.
    1 / (1 - sequencing saturation) can be used can be used to estimate how many additional reads would be required to identify a new UMI. For example, if the sequencing saturation is 0.5 
    (50%) then for each two new reads, one of those should represent a new UMI. If the sequencing saturation is lower, at 0.2
    for example, then on average 1.25 reads would need to be sequenced to obtain a new UMI.


## Make UMAP plots
UMAP (Uniform Manifold Approximation and Projection) is a data visualisation algorithm used
for dimensionality reduction.  It aims to preserve the local structure and relationships in high-dimensional data (in this case a gene x cell count matrix)
when projecting it into a two-dimensional space. This allows structure within the data, such as cell type and state, to be visualised, 
and can be useful for quality control and gaining an initial view of the data.  To generate UMAP plots use `--plot_umaps`. The data used for the 
UMAP generation are the [processed expression matrices](#Make-expression-matrices)


Several UMAP plots are created:
- gene expression plots show the gene expression UMAPs with each cell annotated with the mean gene expression count
- mitochondrial expression plots are the same but annotated with mean mitochondrial expression
- single gene expression UMAPs again show the gene expression UMAPs, but each plot is overlaid with a gene count data from a single gene. 
 These genes of interest can be specified by adding gene names to a file with the path spcified by `--umap_plot_genes path/to/gene+names.csv`
- transcript expression plots shows the transcript expression UMAP plots with each cell annotated with the mean gene transcript count per cell

The UMAP algorithm is stochastic, therefore analysing the same data multiple times, using identical parameters, can lead to visually different projections. 
In order to have some confidence in the observed results, it can be useful to run the projection multiple times.
The number of repeated projections can be set with `--umap_n_repeats` (default 6)

