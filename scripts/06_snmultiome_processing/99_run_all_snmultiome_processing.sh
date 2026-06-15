#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_snmultiome_processing.sh
# Run the snMultiome R analysis workflow from the repository root.
#
# This does not run Cell Ranger ARC by default because that step requires
# large FASTQ files and an HPC/local installation. To run Cell Ranger ARC,
# first run scripts 01 and 02 manually.

Rscript scripts/06_snmultiome_processing/03_load_cellranger_arc_seurat_signac.R
Rscript scripts/06_snmultiome_processing/04_qc_filter_and_merge.R
Rscript scripts/06_snmultiome_processing/05_integrate_wnn_cluster_umap.R
Rscript scripts/06_snmultiome_processing/06_marker_based_celltype_annotation.R
Rscript scripts/06_snmultiome_processing/07_calculate_gene_activity.R
Rscript scripts/06_snmultiome_processing/08_summarize_rna_activity_by_group.R
Rscript scripts/06_snmultiome_processing/09_extract_erf_and_target_summaries.R
Rscript scripts/06_snmultiome_processing/10_qc_check_snmultiome_outputs.R

echo "snMultiome processing workflow complete."
