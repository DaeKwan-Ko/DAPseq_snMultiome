# 02_dapseq_annotation

This folder contains scripts for annotating JGI-provided DAP-seq peak files to the sorghum v3.1.1 genome and for defining conservative DAP-seq-associated candidate target genes.

These scripts start from BED-formatted peak files generated in `scripts/01_dapseq_processing/`, not from raw FASTQ files. JGI-provided DAP-seq peak files are in narrowPeak-style format, where column 7 corresponds to peak summit fold-enrichment. This folder focuses on the downstream annotation workflow used in the manuscript.

## Expected inputs

Edit `00_config.R` before running.

Required inputs:

- Sorghum genome annotation GFF3 file, for example:
  - `Sbicolor_454_v3.1.1.gene.gff3`
- BED files prepared for ChIPseeker, one per TF:
  - expected pattern: `*_for_ChIPseeker.bed`
- Metadata table:
  - `metadata/tf_list_142.tsv`
  - columns: `geneID`, `tf_family`

Optional inputs:

- ERF clade annotation:
  - `metadata/erf_clade_annotation.tsv`
  - columns: `geneID`, `clade`
- Gene-to-GO mapping:
  - `metadata/sorghum_gene2go.tsv`
  - required columns: `geneID`, `go_id`
  - optional columns: `go_name`, `ontology`

## Outputs

Default output directory:

```text
results/02_dapseq_annotation/
```

Main outputs:

- `sorghum_v3.1.1.txdb.sqlite`
- `all_peak_annotations.tsv`
- `peak_annotation_summary_by_tf.tsv`
- `dapseq_candidate_targets_by_tf.tsv`
- `dapseq_candidate_targets_long.tsv`
- `dapseq_target_count_by_tf.tsv`
- `erf_target_matrix.tsv`
- `erf_target_clade_counts.tsv`
- `erf_clade_biased_recurrent_targets.tsv`
- `go_enrichment_by_tf.tsv`
- `go_enrichment_erf_clade_targets.tsv`

## Conservative candidate target definition

DAP-seq-associated candidate targets are defined as genes associated with peaks in one of the following categories:

1. promoter-proximal peaks within 5 kb upstream of a transcription start site;
2. peaks overlapping gene bodies;
3. peaks located within 300 bp downstream of annotated genes.

Distal intergenic peaks are retained in the full annotation table but excluded from GO-based candidate target assignment because their target genes cannot be confidently inferred from genomic proximity alone.

## Suggested run order

```bash
Rscript 01_make_txdb.R
Rscript 02_annotate_peaks_chipseeker.R
Rscript 03_summarize_peak_annotations.R
Rscript 04_define_candidate_targets.R
Rscript 05_go_enrichment_hypergeometric.R
Rscript 06_build_erf_target_matrices.R
```

or run:

```bash
bash 99_run_all_annotation.sh
```

## Notes

Final figure layout was assembled outside this repository. These scripts are intended to reproduce the main peak annotation and candidate target tables used for downstream analyses.

