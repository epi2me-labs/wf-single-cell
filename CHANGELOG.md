# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [v2.0.0]
### Fixed
- Reported cell count off by -1 in report summary table.
- Issue with TSV concat/splitting during `combine_bam_and_tags` stage.
- Issue introduced in v1.1.0 that caused a partial BAM file to be output.
- Corrected example command in README.
- Incorrect reporting of unique gene and transcripts in report table.
- Processed expression matrix entries incorrectly filtered.
- Gene identity of multimapping reads could be incorrectly assigned.
### Changed
- Read chunking done in library code.
- `--process_chunk_size` parameter changed to `--fastq_chunk`
- Resource declarations in Nextflow processes.
- Simplified read batching and decoupled from CPU usage parameters.
- Expression matrix construction code reworked to reduce memory usage.
- Adapter search step now 3x faster.
- Barcode assignment 3x faster.
- Feature assignment now 15x faster.
- UMI clustering 20x faster.
- UMAP creation memory use reduced 6-fold and up-to 30x faster (and
  always enabled).
- Final read tagging step is 3x faster.
- Combined various preprocessing steps into a single process to avoid
  unnecessary file writes.
- Updated stringtie2 to v2.2.2.
- Pre-calculate report summary data to reduce disk-space and IO overheads.
- Single BAM per-sample is now always produced (option `--merge_bam` is removed).
### Removed
- Several workflow parameters as part of resource management simplification.
- `--plot_umaps` option, as UMAP generation has been made much more efficient
  and is always enabled.
- `--merge_bam` option.

## [v1.1.0]
### Added
- `full_length_only` parameter to process only full length reads (default: true).
- Trim adapters, barcodes and UMIs from reads before alignment.
- Memory directive for umap process to prevent parallel processes from using too much memory.
### Changed
- Orient 3prime/multiome reads to mRNA sense to avoid need to flip later.
- Default `umap_n_repeats` lowered to 3.
- Genome reference alignment done by chunk.
### Fixed
- Issue where splice junctions were searched for on incorrect strand.

## [v1.0.3]
### Added
- Publish stringtie transcriptome fasta and GFF files to output dir.
### Fixed
- More informative error message upon read duplicate detection.
### Updated
- Remove duplicate fastcat call.

## [v1.0.2]
### Fixed
- Error interpreting CSV data types during BAM tagging.

## [v1.0.1]
### Fixed
- `<img>` tags in the docs.

## [v1.0.0]
### Updated
- Docs to the new format.

## [v0.3.0]
### Fixed
- `single_cell_sample_sheet` samples with same kit name and version not compatible.

### Changed
-`exp_cells` to `expected_cells` in single_cell_sample_sheet to be consistent with CLI option.

## [v0.2.9]
- Make `prepare_report_data` process more memory-efficient 

## [v0.2.8]
### Fixed
- Increase the maximum memory available to the adapter_scan process
- Fix sequence truncation by 1 bp in adapter_scan step
- Make `summarize_adapter_table` process more memory-efficient

## [v0.2.7]
### Fixed
- Mitochondrial expression file not being copied to output directory
- Incorrect setting of polars maximum threads

### Added
- Allow `geneName` attribute in GTF annotation file

## [v0.2.6]
### Fixed
- Alignments generated from 5' 10x kit are now in the correct orientation.

## [v0.2.5]
### Added
- Memory directives to some processes to better manage system resources 

### Changed
- Bumped minimum required Nextflow version to 22.10.8
- GitHub issue templates
- Add chunking of input data to some processes to reduce memory usage

### Fixed
- Output BAM files with alignments from incorrect chromosomes
- Incorrect uncorrected_barcodes.tsv output

## [v0.2.4]
### Added
- Configuration for running demo data in AWS

## [v0.2.3]
### Fixed
- Barcode assignment error when chromosome has no no data

### Changed
- Include reads in gene expression matrices (but not transcript matrices) that map to intron-only regions

## [v0.2.2]
### Fixed
- Incorrect UMAP colors
- Barcode quality extract error
- Saturation plotting error
- Gene ID assigned instead of gene name
- Empty dataframe bug when no data for a chromosome exists

### Changed
- Improved isoform selection criteria

## [v0.2.1]
### Changed
- Add multiprocessing to calc_saturation.py for ~ 2X speedup
- Use rapidfuzz for finding barcode matches in whielist; ~10x speedup
- Use multithreading to speed up sequencing saturation calcualtion
- Put UMAPs in report and make optional
- ### Changed
- Combine barcode and umi extraction into single step.
- Get gene assignments from stringtie.

## [v0.2.0]
### Added
- workflow-glue to allow scripts to be run as a module.
- pytest testing using workflow container.

### Fixed
- Incorrectly stranded reads causing stringtie2 to generate incorrect transcripts.

## [v0.1.9]
### Fixed
- Incorrect UMIs reported and not collapsing into unique UMI counts.

## [v0.1.8]
### Fixed
- Sample_sheet format in schema to expect a file

## [v0.1.7]
### Changed
- Updated description in manifest.

## [v0.1.6]
### Added
- A workflow report using ezcharts.

### Changed
- Updates for the new Labs Launcher.

## [v0.1.5]
### Changed
- Replace samlmon for minimap2 for assigning reads to transcripts.

### Added
- Reduced matrix to the top N principal components pripor to umap generation.

## [v0.1.4]
### Fixed
- Fix transcript matrices not in output folder.

### Added
- output of merged bam optional.
- Repeat umap creation with different random states.

### Changed
- Transcript counting Salmon on stringtie-generated transcriptome.
- Several performance-related reforactorings including reductions in read write operations.
- single_cell_sample_sheet is optional and kit options can be supplied as workflow parameters.

## [v0.1.3]
### Changed
- Better handling of sample_id conflicts in single_cell_sampkle_sheet and fastgingress.
- single_cell_sample_sheet is optional.
- Minor IO performance enhancements.

### Added
- kit options can be supplied from command line/config and applied to all samples.

## [v0.1.2]
### Added
- Transcript x cell matrix output.

## [v0.1.1]
### Changed
- Combined gathering and splitting of fastqs into a single process.
- Use split2 for splitting fastqs.
- Remove unused kneeplot flags.

### Added
- check for identical sample_ids in single cell sample sheet and fastq data.

## [v0.1.0]
### Added
- First release. Port of [Sockeye](https://github.com/nanoporetech/sockeye) to Nextflow

