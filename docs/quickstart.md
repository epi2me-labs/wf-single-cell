## Quickstart

The workflow uses [nextflow](https://www.nextflow.io/) to manage compute and 
software resources, as such nextflow will need to be installed before attempting
to run the workflow.

The workflow can currently be run using either
[Docker](https://www.docker.com/products/docker-desktop) or
[singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) to provide isolation of
the required software. Both methods are automated out-of-the-box provided
either docker or singularity is installed.

It is not required to clone or download the git repository in order to run the workflow.
For more information on running EPI2ME Labs workflows [visit out website](https://labs.epi2me.io/wfindex).

**Workflow options**

To obtain the workflow, having installed `nextflow`, users can run:

```
nextflow run epi2me-labs/wf-single-cell --help
```

to see the options for the workflow.


The main options are:
* `fastq`: A fastq file or directory containing fastq input files or directories of input files.
* `ref_genome_dir` The path to the 10x reference genome directory (see `Downloading reference data` below)
* 10x sample metadata using either:
  * The following parameters, which are applied to all samples (the default):
    * kit_name 
      * options: `3prime` (default), `5prime`, `multiome`
    * kit_version
      * 3prime options: `v2`, `v3` (default)
      * 5prime options: `v1`
      * multiome options: `v1`
    * expected_cells [500]
* `plot_umaps`: This flag controls whether UMAP projections are generated and dispalyed in the report (default `false`). If UMAP output is required apply like: `--plot_umaps`. 

  or

  *  `single_cell_sample_sheet`
  (__not to be confused with the optional MinKNOW `sample_sheet`__) allowing per sample configuration.



The single_cell_sample_sheet contains details about the input sample_ids, the 10X kits used (e.g. `3prime` or `5prime`), the kit versions used (`v2` or `v3` for the 3' kit, `v1` for the 5' kit), a rough estimate of the number of cells in the library. The cell count estimate specified with `exp_cells` and can be a very rough estimate (500 is a robust default value if the number is not known).


The sample_id field should correspond to sample_id which is defined either in the `sample_sheet`,  given by the `sample` parameter (for single sample runs). If no `sample_sheet` or `sample` is given, sample_id is derived from each folder containing the fastq files or if a single file is given, the sample_id is the basename of the file (data.fastq.gz -> data.fastq).

An example sheet with one sample is:
```
sample_id,kit_name,kit_version,exp_cells
sample_10,3prime,v3,500
```

**Downloading reference data**
The pipeline requires access to reference data files that are packaged and freely available from 10x Genomics. For human samples, the GRCh38 packaged reference files can be downloaded using either curl or wget using:

```
curl -O https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

or 
```
wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2020-A.tar.gz
tar -xvf refdata-gex-GRCh38-2020-A.tar.gz
```

**Download demonstration data**

A dataset of around 500 MB is provided for the purposes of testing the workflow.
It consists of a downsampled dataset of reads mapping to human chromosome 17. 
The biological sources of the data are five cancer cell lines ([see this article](https://www.nature.com/articles/s41592-019-0425-8))
, which can be visualied as clusters in the resulting UMAP plots.
It can be downloaded using:

```
wget https://ont-exd-int-s3-euwst1-epi2me-labs.s3.amazonaws.com/wf-single-cell/wf-single-cell-demo.tar.gz \
  && tar -xzvf wf-single-cell-demo.tar.gz
```

The workflow can be run with the demonstration data using:

```
OUTPUT=output
nextflow run epi2me-labs/wf-single-cell \
    -w ${OUTPUT}/workspace \
    -profile standard \
    --fastq wf-single-cell-demo/chr17.fq.gz \
    --kit_name 3prime \
    --kit_version v3 \
    --expected_cells 100 \
    --ref_genome_dir wf-single-cell-demo/ \
    --out_dir ${OUTPUT} \
    --plot_umaps
```

The output of the pipeline will be found in `./output` for the above
example. 

**Supplying your own reference data**

To use your own reference data, a folder should be created following the structure
of `reference_data` shown below.

```
reference_data
├── fasta
│   ├── genome.fa
│   └── genome.fa.fai
└── genes
    └── genes.gtf
```
Note that the annotation file should be in GTF format and the transcript feature lines should contain the attributes `gene_name` or `geneName` as well as the attribute `transcript_id`.

The reference data folder can then be referred to by the `--ref_genome_dir` parameter.

**Workflow outputs**

The pipeline output will be written to a directory defined by ``--out_dir``. 
For each sample specified in the `single_cell-sample_sheet`  an output folder is created containing the results.


The most useful outputs of the pipeline are likely:

* ``configs.stats.json``: provides a summary of sequencing statistics and observed read configurations, such as

  - ``n_reads``: number of total reads in the input fastq(s)
  - ``rl_mean``: mean read length
  - ``n_fl``: total number of reads with the read1-->TSO or TSO'-->read1' adapter configuration (i.e. full-length reads)
  - ``n_plus``: number of reads with the read1-->TSO configuration
  - ``n_minus``: number of reads with the TSO'-->read1' configuration

* ``bams``: Folder of bam alignment files where each alignment contains the following sequence tags

  - CB: corrected cell barcode sequence
  - CR: uncorrected cell barcode sequence
  - CY: Phred quality scores of the uncorrected cell barcode sequence
  - UB: corrected UMI sequence
  - UR: uncorrected UMI sequence
  - UY: Phred quality scores of the uncorrected UMI sequence


  The bam files are output per chromosome (default) unless `--merge_bam` is set.

 * ``gene_expression.processed.tsv``:  TSV containing the gene (rows) x cell (columns) expression matrix, processed and normalized according to: 

  - ``matrix_min_genes``: cells with fewer than this number of expressed genes will be removed
  - ``matrix_min_cells``: genes present in fewer than this number of cells will be removed
  - ``matrix_max_mito``: cells with more than this percentage of counts belonging to mitochondrial genes will be removed
  - ``matrix_norm_count``: normalize all cells to this number of total counts per cell


* ``transcript_matrix_processed.tsv``: TSV containing the transcript (rows) x cell (columns) expression matrix, processed and normalized in the same manner as the genes.
These expression values are determined by first generating a transcriptome per sample using [stringtie](https://ccb.jhu.edu/software/stringtie/) and then assigning reads to
transcripts aligning them to this transcriptome with minimap2. Only reads that map unambiguously to a reference transcript are assigned. 
The assembled transcripts with the following gffcompare class codes
are excluded: `i`, `p`, `s` or `u`.
See the [gffcompare](https://ccb.jhu.edu/software/stringtie/gffcompare.shtml)
and [this image](https://ccb.jhu.edu/software/stringtie/gffcompare_codes.png), and
only cells and genes that pass the gene filtering described above are included.
* ``read_tags.tsv``: 
TSV file witjh the following columns:
  - read_id
  - gene (assigned gene)
  - transcript (assigned transcript id)  
  - barcode (corrected barcode)
  - umi )corrected umi