# 06_snmultiome_processing

This folder contains scripts for processing and analyzing sorghum root single-nucleus multiome (snMultiome) data.

The workflow documents Cell Ranger ARC processing, Seurat/Signac quality control, RNA–ATAC integration, marker-based cell-type annotation, gene activity calculation and ERF expression/activity summaries used for cell type-resolved salinity analyses.

## Main goals

1. Document Cell Ranger ARC reference building and library processing.
2. Load control and NaCl-treated 10x Multiome outputs into Seurat/Signac.
3. Apply joint RNA and ATAC quality-control filtering.
4. Integrate RNA and ATAC modalities using weighted nearest neighbors.
5. Annotate major sorghum root cell populations using marker modules.
6. Calculate RNA expression and ATAC-derived gene activity summaries by cell type and condition.
7. Summarize expression and gene activity of DAP-seq-profiled ERFs and candidate target genes.

## Expected metadata

Required:

```text
metadata/sample_metadata_snmultiome.tsv
```

Recommended columns:

```text
sample_id
condition
treatment
treatment_duration_h
genotype
species
tissue
assay
modalities
cellranger_arc_version
reference_genome
reference_annotation
bioproject
biosample
sra_atac_accession
sra_gex_accession
nuclei_retained_after_qc
```

Also used by downstream scripts when available:

```text
metadata/erf_clade_annotation.tsv
metadata/salinity_linked_erfs.tsv
metadata/tf_list_142.tsv
metadata/root_marker_genes.tsv
metadata/candidate_target_genes.tsv
```

`root_marker_genes.tsv` expected columns:

```text
celltype
geneID
marker_source
```

`candidate_target_genes.tsv` expected columns:

```text
geneID
target_group
```

## Expected Cell Ranger ARC outputs

Place Cell Ranger ARC output folders here or adjust `00_config.R`:

```text
results/06_snmultiome_processing/input/cellranger_arc/
├── Ctl1/outs/
└── NaCl1/outs/
```

Each `outs/` directory should contain standard Cell Ranger ARC files such as:

```text
filtered_feature_bc_matrix.h5
atac_fragments.tsv.gz
per_barcode_metrics.csv
```

Large Cell Ranger output files should not be committed to GitHub. These paths document where files are expected after downloading raw data from SRA and rerunning Cell Ranger ARC locally or on HPCC.

## Default QC thresholds

The default thresholds match the manuscript workflow:

```text
RNA counts: 200–50,000
ATAC counts: 200–50,000
TSS enrichment: >1.2
Nucleosome signal: <0.35
```

## Outputs

Default output directory:

```text
results/06_snmultiome_processing/
```

Main outputs:

- `snmultiome_merged_qc_filtered.rds`
- `snmultiome_wnn_integrated.rds`
- `celltype_annotation_metadata.tsv`
- `celltype_condition_nucleus_counts.tsv`
- `gene_activity_matrix.rds`
- `rna_by_celltype_condition.tsv`
- `activity_by_celltype_condition.tsv`
- `erf_rna_activity_summary_by_celltype_condition.tsv`
- `target_gene_rna_activity_summary_by_celltype_condition.tsv`
- `snmultiome_processing_qc.tsv`

## Suggested run order

From the repository root:

```bash
bash scripts/06_snmultiome_processing/01_build_cellranger_arc_reference.sh
bash scripts/06_snmultiome_processing/02_run_cellranger_arc_count.sh
Rscript scripts/06_snmultiome_processing/03_load_cellranger_arc_seurat_signac.R
Rscript scripts/06_snmultiome_processing/04_qc_filter_and_merge.R
Rscript scripts/06_snmultiome_processing/05_integrate_wnn_cluster_umap.R
Rscript scripts/06_snmultiome_processing/06_marker_based_celltype_annotation.R
Rscript scripts/06_snmultiome_processing/07_calculate_gene_activity.R
Rscript scripts/06_snmultiome_processing/08_summarize_rna_activity_by_group.R
Rscript scripts/06_snmultiome_processing/09_extract_erf_and_target_summaries.R
Rscript scripts/06_snmultiome_processing/10_qc_check_snmultiome_outputs.R
```

or run:

```bash
bash scripts/06_snmultiome_processing/99_run_all_snmultiome_processing.sh
```

## Notes

The scripts are written to document and reproduce the computational workflow. Large raw FASTQ files, Cell Ranger ARC output folders, fragment files and Seurat objects should not be committed to GitHub. Store those in public archives or local/HPC storage and provide accession numbers in the manuscript and metadata.

