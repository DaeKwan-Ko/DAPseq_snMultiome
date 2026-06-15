#!/usr/bin/env bash
set -euo pipefail

# 99_run_all_erf_phylogeny_motifs.sh
# Run the ERF phylogeny and motif workflow from the repository root.
#
# Note:
# The beta2-signature script requires manual alignment-column specification.
# Set BETA2_COLUMNS before running, for example:
#   export BETA2_COLUMNS=25,26
#
# Usage:
#   bash scripts/03_erf_phylogeny_motifs/99_run_all_erf_phylogeny_motifs.sh

source scripts/03_erf_phylogeny_motifs/00_config.sh

Rscript scripts/03_erf_phylogeny_motifs/01_prepare_erf_metadata.R

bash scripts/03_erf_phylogeny_motifs/01_run_hmmer_ap2_scan.sh

python scripts/03_erf_phylogeny_motifs/02_extract_ap2_domains_from_hmmer.py

Rscript scripts/03_erf_phylogeny_motifs/03_align_dbd_and_build_tree.R

if [[ -n "${BETA2_COLUMNS:-}" ]]; then
  python scripts/03_erf_phylogeny_motifs/04_assign_beta2_signature_from_alignment.py \
    --columns "${BETA2_COLUMNS}"
else
  echo "Skipping beta2 signature extraction because BETA2_COLUMNS is not set."
  echo "After checking the alignment, rerun with: export BETA2_COLUMNS=col1,col2"
fi

Rscript scripts/03_erf_phylogeny_motifs/05_parse_meme_top_motifs.R
Rscript scripts/03_erf_phylogeny_motifs/06_integrate_phylogeny_motif_metadata.R
bash scripts/03_erf_phylogeny_motifs/07_make_ceqlogo_commands.sh
Rscript scripts/03_erf_phylogeny_motifs/08_qc_check_erf_outputs.R

echo "ERF phylogeny and motif workflow complete."
