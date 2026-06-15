#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_annotation.sh
# Run the DAP-seq annotation workflow from the repository root.
#
# Usage:
#   bash scripts/02_dapseq_annotation/99_run_all_annotation.sh

Rscript scripts/02_dapseq_annotation/01_make_txdb.R
Rscript scripts/02_dapseq_annotation/02_annotate_peaks_chipseeker.R
Rscript scripts/02_dapseq_annotation/03_summarize_peak_annotations.R
Rscript scripts/02_dapseq_annotation/04_define_candidate_targets.R
Rscript scripts/02_dapseq_annotation/05_go_enrichment_hypergeometric.R
Rscript scripts/02_dapseq_annotation/06_build_erf_target_matrices.R
Rscript scripts/02_dapseq_annotation/07_qc_check_annotation_outputs.R

echo "DAP-seq annotation workflow complete."
