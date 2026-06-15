#!/usr/bin/env bash
set -euo pipefail

# 01_run_hmmer_ap2_scan.sh
# Run HMMER scanning to identify AP2/ERF DNA-binding domains.
#
# Usage from repository root:
#   bash scripts/03_erf_phylogeny_motifs/01_run_hmmer_ap2_scan.sh
#
# Required:
#   PROTEIN_FASTA: protein FASTA containing sorghum ERF or TF sequences
#   HMM_DATABASE:  HMM database containing AP2 domain model, e.g. Pfam-A.hmm

source scripts/03_erf_phylogeny_motifs/00_config.sh

if [[ ! -s "${PROTEIN_FASTA}" ]]; then
  echo "ERROR: PROTEIN_FASTA not found or empty: ${PROTEIN_FASTA}" >&2
  exit 1
fi

if [[ ! -s "${HMM_DATABASE}" ]]; then
  echo "ERROR: HMM_DATABASE not found or empty: ${HMM_DATABASE}" >&2
  exit 1
fi

DOMTBLOUT="${OUT_DIR}/sorghum_erf_hmmer.domtblout"
STDOUT_LOG="${LOG_DIR}/hmmscan_ap2.stdout.txt"

echo "Running HMMER scan..."
echo "Protein FASTA: ${PROTEIN_FASTA}"
echo "HMM database: ${HMM_DATABASE}"

"${HMMSCAN_BIN}" \
  --cpu 8 \
  --domtblout "${DOMTBLOUT}" \
  "${HMM_DATABASE}" \
  "${PROTEIN_FASTA}" \
  > "${STDOUT_LOG}"

echo "Saved domtblout: ${DOMTBLOUT}"
echo "Saved stdout log: ${STDOUT_LOG}"
