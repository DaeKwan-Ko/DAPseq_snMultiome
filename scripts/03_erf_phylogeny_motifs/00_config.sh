#!/usr/bin/env bash

# 00_config.sh
# Shell configuration for ERF phylogeny and motif scripts.
# Edit these paths before running the shell/Python scripts.

# Repository root. The default assumes scripts are run from the repository root.
export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

# Metadata
export METADATA_DIR="${PROJECT_DIR}/metadata"
export TF_METADATA="${METADATA_DIR}/tf_list_142.tsv"
export ERF_CLADE_METADATA="${METADATA_DIR}/erf_clade_annotation.tsv"
export SALINITY_ERF_METADATA="${METADATA_DIR}/salinity_linked_erfs.tsv"

# Input sequence files
# Use a FASTA containing the 49 sorghum ERF protein sequences, or all 142 TF proteins.
export PROTEIN_FASTA="${PROJECT_DIR}/reference/sorghum_dapseq_tf_proteins.fa"

# HMMER database/profile.
# Example: Pfam-A.hmm containing the AP2 domain model PF00847.
export HMM_DATABASE="${PROJECT_DIR}/reference/Pfam-A.hmm"
export HMM_DOMAIN_REGEX="AP2|PF00847|AP2/ERF"

# MEME output directory from DAP-seq motif discovery.
# The scripts recursively search this directory for meme.txt and .meme files.
export MEME_DIR="${PROJECT_DIR}/results/01_dapseq_processing/meme"

# Optional DAP-seq peak summary table.
export PEAK_SUMMARY="${PROJECT_DIR}/results/01_dapseq_processing/peak_count_summary.tsv"

# Output directories
export OUT_DIR="${PROJECT_DIR}/results/03_erf_phylogeny_motifs"
export INTERMEDIATE_DIR="${OUT_DIR}/intermediate"
export LOG_DIR="${OUT_DIR}/logs"
export LOGO_DIR="${OUT_DIR}/motif_logos"

# Executables
export HMMSCAN_BIN="${HMMSCAN_BIN:-hmmscan}"
export CLUSTALW_BIN="${CLUSTALW_BIN:-clustalw2}"
export CEQLOGO_BIN="${CEQLOGO_BIN:-ceqlogo}"

mkdir -p "${OUT_DIR}" "${INTERMEDIATE_DIR}" "${LOG_DIR}" "${LOGO_DIR}"

echo "PROJECT_DIR=${PROJECT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
