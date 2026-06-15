#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_evolutionary_analysis.sh
# Run the evolutionary analysis workflow from the repository root.
#
# Usage:
#   bash scripts/04_evolutionary_analysis/99_run_all_evolutionary_analysis.sh
#
# Optional:
# To generate cross-species beta2 signatures from an aligned FASTA, first run:
#   python scripts/04_evolutionary_analysis/04_assign_beta2_signature_cross_species.py --columns col1,col2
#
# where col1,col2 are manually verified beta2 alignment columns.

source scripts/04_evolutionary_analysis/00_config.sh

Rscript scripts/04_evolutionary_analysis/01_prepare_species_metadata.R
Rscript scripts/04_evolutionary_analysis/02_parse_orthofinder_orthologues.R
Rscript scripts/04_evolutionary_analysis/03_prepare_kaks_ortholog_pairs.R
Rscript scripts/04_evolutionary_analysis/04_summarize_kaks_by_clade.R
Rscript scripts/04_evolutionary_analysis/05_build_cross_species_beta2_summary.R
Rscript scripts/04_evolutionary_analysis/06_prepare_combined_erf_phylogeny_metadata.R
Rscript scripts/04_evolutionary_analysis/07_qc_check_evolutionary_outputs.R

echo "Evolutionary analysis workflow complete."
