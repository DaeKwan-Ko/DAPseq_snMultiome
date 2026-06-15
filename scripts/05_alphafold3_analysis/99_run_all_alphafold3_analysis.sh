#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_alphafold3_analysis.sh
# Run the AlphaFold3 analysis workflow from the repository root.
#
# Usage:
#   bash scripts/05_alphafold3_analysis/99_run_all_alphafold3_analysis.sh
#
# Before running, create:
#   metadata/af3_designs.tsv
#   metadata/af3_motif_position_map.tsv
#   metadata/af3_beta2_residue_map.tsv
#
# and place AlphaFold3 downloaded output folders in:
#   results/05_alphafold3_analysis/input/af3_outputs/

python scripts/05_alphafold3_analysis/01_prepare_af3_inputs.py
python scripts/05_alphafold3_analysis/02_collect_af3_outputs.py
python scripts/05_alphafold3_analysis/03_extract_structure_atom_tables.py
python scripts/05_alphafold3_analysis/04_calculate_dna_rmsd.py
python scripts/05_alphafold3_analysis/05_calculate_protein_dna_contacts.py
Rscript scripts/05_alphafold3_analysis/06_summarize_af3_metrics.R
Rscript scripts/05_alphafold3_analysis/07_qc_check_af3_outputs.R

echo "AlphaFold3 analysis workflow complete."
