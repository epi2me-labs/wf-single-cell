<!---High level numbered list of main steps of the workflow and hyperlink to any tools used. If multiple workflows/different modes perhaps have subheadings and numbered steps. Use nested numbering or bullets where required.--->

The following section details the main steps of the workflow. 

### 1. Concatenates input files and generate per read stats.

The [fastcat/bamstats](https://github.com/epi2me-labs/fastcat) tool is used to concatenate multi-file samples to be processed by the workflow. It will also output per read stats including average read lengths and qualities.

### 2. Stranding and identification of full-length reads

Reads derived from a 10x Genomics library will ideally be flanked by two different adapter
sequences. These reads are more likely to represent full length mRNA molecules,
although that isn't guaranteed as some cDNAs may have been created by internal
priming or represent other cDNA synthesis artifacts. See this [10x Genomics Technical Note](https://cdn.10xgenomics.com/image/upload/v1660261286/support-documents/CG000376_TechNote_Antisense_Intronic_Reads_SingleCellGeneExpression_RevA.pdf) for more information.


This stage of the workflow assigns an adapter configuration to each read. 
This assigned configurations allows the stranding and orientating the reads. 
The following schematic shows an example read structure from 10x Genomics 3&#8242; library .

<figure>
<img src="docs/images/3prime_read.png" alt="10x read structure"/>
<figcaption>Fig.1 Read structure for 10x 3prime kit reads</figcaption>
</figure>

Adapters are located within the reads using [vsearch](https://github.com/torognes/vsearch) (`Read1` and `TSO` in Fig.1 in the case of the 3prime kit).
The table below details the adapter sequences for each of the 10x Genomics kits, along with links to the relevant user guides.


| Kit      | adapter1                | adapter2                     | 10x user guide                                                                                                                              |
|----------|-------------------------|------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------|
| 3&#8242; | CTACACGACGCTCTTCCGATCT  | ATGTACTCTGCGTTGATACCACTGCTT  | [3' kit](https://cdn.10xgenomics.com/image/upload/v1660261285/support-documents/CG000204_ChromiumNextGEMSingleCell3_v3.1_Rev_D.pdf)         |
| multiome | CTACACGACGCTCTTCCGATCT  | ATGTACTCTGCGTTGATACCACTGCTT  | [5'  kit](https://cdn.10xgenomics.com/image/upload/v1666737555/support-documents/CG000331_ChromiumNextGEMSingleCell5-v2_UserGuide_RevE.pdf) |
| 5&#8242; | CTACACGACGCTCTTCCGATCT  | GTACTCTGCGTTGATACCACTGCTT    | [Multiome kit](https://teichlab.github.io/scg_lib_structs/data/CG000338_ChromiumNextGEM_Multiome_ATAC_GEX_User_Guide_RevB.pdf)              |
| Visium   | CTACACGACGCTCTTCCGATCT  | ATGTACTCTGCGTTGATACCACTGCTT  | [Visium kit](https://cdn.10xgenomics.com/image/upload/v1695417753/support-documents/CG000239_VisiumSpatialGeneExpression_UserGuide_RevG.pdf) |

Concatenated reads are identified by adapter configuration and are split into individual subreads and 
reorientated if required.
The following table details the various configurations and the actions taken for each.


| configuration   | action                                                    |
|-----------------|-----------------------------------------------------------|
| Full length     | Reads are trimmed from each side and oriented and adapter2 is removed           |
| Single adapter1 | Reads are oriented and trimmed from adapter1 end only     |
| Single adapter2 | Reads are oriented and trimmed from the adapter2 end only |
| double adapter1 | Reads are trimmed from both sides                         |
| double adapter2 | Reads are trimmed from both sides                         |
| other           | No valid adapters found; not used in further analysis     | 

Adapter configuration summaries can be found in the output file  `{{ alias }}/{{ alias }}.config_stats.json"`

To only process full length reads the option `--full_length_only` should be set to true (default: true). 
If set to false, reads with only a single adapter or other non-full-length adapter configurations will also be processed.


### 3. Extract cell barcodes and UMIs
The next step is to extract 10x Genomics barcodes and UMI sequences from the stranded and trimmed reads.

In order to do this, the first 100bp of each read are aligned to a reference probe using [parasail](https://github.com/jeffdaily/parasail). This probe contains a suffix of the adapter1 sequence, some ambiguities ("Ns") representing the barcode and UMI, and a polyT tract.

<figure>
<img src="docs/images/probe.png" alt="probe image"/>
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
- The size of the adapter1 suffix can be specified with: `barcode_adapter1_suff_length`

Once the barcode and UMI has been extracted, these are then trimmed from the reads along
with the adapter1 sequence. These trimmed reads are then used in the alignment step.

### 4. Aligning reads to genome
The next stage is to align the preprocessed reads to the reference genome. This enables gene and transcript
read assignment in downstream steps.

The stranded and trimmed FASTQ reads are mapped to the reference genome using minimap2.  

The reference data can be supplied one of two ways:
1. as a local path to a folder with `--ref_genome_dir`.
Either use a 10x reference bundle (https://www.10xgenomics.com/support/software/cell-ranger/downloads#References) 
or ensure that the reference directory contains the following files in the same folder structure
```
├── refdata
│   ├── fasta
│   │   ├── genome.fa
│   │   └── genome.fa.fai
│   └── genes
│       └── genes.gtf
```
2. `--epi2me_resource_bundle`: Select this option to use a prebuilt 10x resource directory.
This will be downloaded automatically by the workflow and stored in `store-dir`;
subsequent runs from the same directory will reuse the reference data stored there.  
There are currently prebuilt resource bundles for the 10x Human reference (GRCh38) - 2024-A
reference bundle (https://www.10xgenomics.com/support/software/cell-ranger/downloads).




### 5. Barcode correction
The aim of this stage is to correct errors present in the previously extracted barcodes.
10x Genomics barcodes are not random; all possible barcode sequences can be found in a whitelist of known barcodes. 
The appropriate whitelist for each kit and version is chosen automatically by the workflow.
The whitelist is used to generate a shortlist of high quality barcodes present in the sample, which is then used to 
correct barcode errors.

The correction proceeds as follows:
* A shortlist of all high quality barcodes present in our sample is generated. 
This is done by adding an uncorrected barcode to the shortlist if it:
  * has a min quality > `barcode_min_quality` (default 15)

In each cell library there are expected to be some low quality cells and empty droplets that can be identified by their low number of reads.
To remove these cells, the shortlist is filtered with a quantile-based threshold if `estimate_cell_count` is set to `true` (default).
This threshold is determined by ranking the cells by read count and taking the top n cells (n = `expected_cells`).
The read count 95th percentile / 20 is the threshold used. This threshold can be visualised in the knee plots generated by the workflow.

Alternatively if a predetermined number of cells is required for analysis, setting `estimate_cell_count` to `false` results in a cell count of `expected_cells`.

*note*: As `visium` barcodes do not represent cells, but rather tissue coordinates,
shortlist cell count thresholding is not performed for `visium` analysis.

For uncorrected barcodes not present in the shortlist, they are cross-referenced against the shortlist, and are assigned 
a barcode from this list if the candidate barcode meets the following criteria:
* the query and closest-matched shortlist barcode have an edit distance <= 2
* The edit distance difference between the top and hit and second top hit in the shortlist >= 2

*note*: Edit distance here refers specifically to the Levenshtein distance, which is the minimum number of single base 
changes it would take to transform one UMI into another (including insertions, substitutions or deletions).

Workflow options:
* `barcode_max_ed`: Max edit distance between the uncorrected barcode and the matching whitelist barcode (default 2).
* `barcode_min_ed_diff`: Min difference in edit distance between (1) the uncorrected barcode vs top hit and (2) uncorrected barcode vs runner-up hit (default 2).
* `expected_cells`: Number of expected cells. Enter the number of expected cells in the sample.


### 6. Gene and transcript assignment
Now that barcodes and UMIs have been assigned, the next stage is to assign a gene and transcript to each read. 

The first step is to create a transcriptome using stringtie2 (using long read mode `-L`).
The minimum required transcripts for a transcript to be called is 2. This can be changed, along with other stringtie settings using the 
`stringtie_opts` options (default "-c 2"). The resulting query transcriptome is then annotated with the supplied reference 
genome annotation using [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml), annotating query transcripts with reference transcript and gene metadata.

The original reads are then mapped to this transcriptome, producing one or more alignments per read.

#### Transcript assignment
To assign a read to a transcript, we follow a similar procedure to that used by [FLAMES](https://github.com/LuyiTian/FLAMES).
* Transcripts that map to intronic query transcripts (those mapping to query transcripts with gffcompare classes ['i', 'y', 'p', 's'])
are set to unknown (-). This does not affect the gene assignment.
* Any read that uniquely maps to the annotated transcriptome is assigned that transcript.
* Any read that does not uniquely map to a transcript is processed with the following procedure:
  * Order alignments by alignment score (AS)
  * If the top alignment has a better AS or query coverage than the second, and
  has a minimum transcript coverage of 0.4, assign to this transcript, else set as unknown (-)
  * If the AS is the same between the 1st and 2nd alignment the choose the alignment that has the highest transcript coverage, if the coverage is > 0.8 

#### Gene assignment
For each read, the alignment with the highest genomic MAPQ score is used for gene assignemnt if this score > `gene_assigns_minqv` (default 30).
The gene ID is then derived from the transcript assigned in the transcript assignment step.

### 7. UMI correction 
The next step is to correct UMI errors. The previous barcode correction step leveraged a
whitelist of all potential barcodes, drastically narrowing the search space. However, UMIs are random sequences and so there is no way of knowing beforehand
which UMIs will be in our library. To reduce the UMI search space, the reads are initially clustered by assigned gene, or if
no gene is assigned, by genomic interval. This also has the effect of reducing the likelihood of UMI collisions which occur due to the number
of total reads in a library usually being significantly higher than the total possible number of UMIS (16.7 million combinations for a 12bp UMI).

UMI clusters are then generated from these gene-clustered reads using [UMI-tools](https://github.com/CGATOxford/UMI-tools) with 
a slightly modified version of the [directional method](https://umi-tools.readthedocs.io/en/latest/the_methods.html).
For this clustering, Levenshtein distance threshold is used (default 2) instead of Hamming distance, in order to account for UMI indels. 
The selected node of each cluster is the one with the highest number of reads and is used to assign a `corrected_umi` to the rest of the reads in that cluster.

### 8. Make expression matrices
The gene x cell and transcript x cell expression matrices are the main outputs of the workflow and can be used in further analysis of 
the single cell experiment. 

The expression matrices are generated by collapsing the corrected UMIs (counting each unique UMI once) and summing
the counts of features (gene or transcript) per cell to give a feature x cell expression matrix and are output as a folder of files in  Market Exchange (MEX) 
format (see the [output docs](docs/outputs.md)  

The expression count matrices are further processed in the following way to give gene x cell processed matrices 
and are also output in  Market Exchange (MEX) format.
* Cells are dropped that contain less than `matrix_min_genes` genes or transcripts (default 200)
* Genes are dropped which are present in fewer than `matrix_min_cells` (default 3)
* Cells where mitochondrial genes make up more than `matrix_max_mito` (default 20%) are dropped
* Counts are normalized to `matrix_norm_count` (default 10,000) reads/cell
* Normalized counts are finally log10 transformed


### 9. Tagging bam files
BAM files generated from aligning reads to the reference genome are now tagged with the 
following information from the workflow. By default, the BAMs are output per chromosome, but can be concatenated
into a single file per sampe using `merge_bam`. The new bam will contain the following tags:
  - CB: corrected cell barcode sequence
  - CR: uncorrected cell barcode sequence
  - CY: Phred quality scores of the uncorrected cell barcode sequence
  - UB: corrected UMI sequence
  - UR: uncorrected UMI sequence
  - UY: Phred quality scores of the uncorrected UMI sequence


### 10. Calculate library saturation
Sequencing saturation is an estimate of the library complexity that has been captured in a sequencing run.
As read depth per cell increases, the number of genes and distinct UMIs identified will increase at a rate that is 
dependent on the complexity of the input library. 

Reads which have been assigned a corrected barcode, corrected UMI and gene
are subsampled by varying degrees and the resulting median reads/cell is plotted against either:
- median genes per cell: this gives an indication of the gene complexity captured and whether increasing read depth would lead 
 to significantly more genes being identified.
- median UMIs per cell: this gives an indication of the UMI complexity captured. If this plot plateaus out, it may indicate that 
    PCR duplicates represent a high proportion of the reads at higher sequencing depth.
- Sequencing saturation:  calculated as `1 - (number of unique UMIs across all cells / number of reads across all cells)`. Values near 0 indicate that
 very few duplicate UMIs are being identified, whereas values nearer to 1 indicate a higher proportion of duplicate UMIs are being captured.
    1 / (1 - sequencing saturation) can be used can be used to estimate how many additional reads would be required to identify a new UMI. For example, if the sequencing saturation is 0.5 
    (50%) then for each two new reads, one of those should represent a new UMI. If the sequencing saturation is lower, at 0.2
    for example, then on average 1.25 reads would need to be sequenced to obtain a new UMI.


### 11. SNV calling
wf-single-cell contains an experimental single nucleotide variant (SNV) calling workflow based on [longshot](https://github.com/pjedge/longshot). 
Currently the workflow has been tested with a maximum of 1500 cells, and can be expected to take up 24 hours with a 64 core machine. 
Work is underway to improve the performance of the SNV workflow.

#### Summary of the SNV workflow:
- Tagged BAMs are split by cell barcode.
- Preprocessing of the barcode-split BAM is required before processing with longshot:
  * UMI deduplication with [UMI-tools](https://umi-tools.readthedocs.io/en/latest/reference/dedup.html).
  * Splitting of records by exon using [gatk SplitNCigarReads](https://gatk.broadinstitute.org/hc/en-us/articles/360036858811-SplitNCigarReads).
- Each per-cell BAM is processed with longshot to generate an initial set of per-cell candidates.
These candidates can represent variants that may not be detected in bulk cDNA as they may be present in few cells.
- The per-cell BAMs for each sample are merged and processed with longshot to give a set of bulk sample candidates. 
Some of these candidate variants may represent variants that were not detected at the initial per-cell genotyping due to low coverage within the cell.
- The bulk sample and per-cell candidate variants are merged to produce a final set of candidate variants.
- A second round of per cell genotyping is carried out with longshot using the merged candidates from the previous steps. 
- The SNV workflow outputs: 
  * A merged VCF (`output/<alias>.final_merged.vcf.gz`) containing the genotype calls for each cell in the sample columns.
  * A MEX format genotype snv x cell matrix (`output/genotype_matrix`), which can be loaded into downstream tools such as seurat.
  * The genotypes are encodes as follows
      * homozygous REF : 0
      * heterozygous ALT/REF 1
      * homozygous ALT: 2


### 12. Fusion transcript calling
Fusion transcript calling can be enabled using [ctat-LR-fusion](https://github.com/TrinityCTAT/CTAT-LR-fusion),
allowing reads derived from fusion transcripts to be identified and assigned to cells. 
This part of the workflow can be enabled with `--call_fusions`.

The taqged BAM output from the workflow, containing cell and UMI barcode tags, is used as
input to ctat-LR-fusion. 

ctat-LR-fusion requires a resource directory containing reference sequence and annotation information. It is important that the ctat-LR-fusion resource is built against the same
reference data as is used elsewhere in the workflow. For example, if the `ref_genome_dir`
contains sequence from hg38 and Gencode44 annotations, then the ctat-LR-fusion resource directory should be built against these.

There are two ways to supply the ctat-LR-fusion resource directory:

1.  `--ctat_resource_dir`: A path to a local copy of the resource directory.
Prebuilt resource directories can be found here:  https://data.broadinstitute.org/Trinity/CTAT_RESOURCE_LIB/.
2. `--epi2me_resource_bundle`: Select this option to use a prebuilt ctat-LR-fusion resource bundle and the corresponding 10x reference data,
which will be automatically downloaded from the cloud.
Currently we have prebuilt bundles for the 10x human reference (GRCh38) - 2024-A
reference bundle (https://www.10xgenomics.com/support/software/cell-ranger/downloads).

The main output is a per read summary file where each read called as a fusion
by ctat-LR-fusion is associated with a cell barcode/UMI and gene/transcript assignments. 


### 12. Make UMAP plots
UMAP (Uniform Manifold Approximation and Projection) is a data visualisation algorithm used
for dimensionality reduction. It aims to preserve the local structure and relationships in high-dimensional data (in this case a gene x cell count matrix)
when projecting it into a two-dimensional space. This allows structure within the data, such as cell type and state, to be visualised, 
and can be useful for quality control and gaining an initial view of the data.  To generate UMAP plots use `plot_umaps`. The data used for the 
UMAP generation are the [processed expression matrices](#Make-expression-matrices)


Several UMAP plots are created:
- gene expression plots show the gene expression UMAPs with each cell annotated with the mean gene expression count
- mitochondrial expression plots are the same but annotated with mean mitochondrial expression
- single gene expression UMAPs again show the gene expression UMAPs, but each plot is overlaid with a gene count data from a single gene. 
 These genes of interest can be specified by adding gene names to a file with the path spcified by `umap_plot_genes path/to/gene+names.csv`
- transcript expression plots shows the transcript expression UMAP plots with each cell annotated with the mean gene transcript count per cell

The UMAP algorithm is stochastic, therefore analysing the same data multiple times, using identical parameters, can lead to visually different projections. 
In order to have some confidence in the observed results, it can be useful to run the projection multiple times.
The number of repeated projections can be set with `umap_n_repeats` (default 3)


### 13 Visium HD support
Visium HD is supported. ONT reads must first be processed by 10x Genomics Space Ranger. 
Please see the instructions at https://epi2me.nanoporetech.com/epi2me-docs/tools/percula/. 

ONT long read BAMs are processed by Percula to produce BAM files that are compatible with Space Ranger.

wf-single-cell has the following relevant options:

- `--bam`: the long read BAM output by Percula
- `--spaceranger_bam`: the demultiplexed tagged BAM output by Space Ranger
- `--adapter_stats`: the configs.json file output by Percula

Note that the workflow expects a single sample for analysis of Visium HD data.