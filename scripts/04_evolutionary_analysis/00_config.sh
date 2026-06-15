#!/usr/bin/env bash

# 00_config.sh
# Shell configuration for evolutionary analysis helper scripts.

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

export METADATA_DIR="${PROJECT_DIR}/metadata"
export ERF_CLADE_METADATA="${METADATA_DIR}/erf_clade_annotation.tsv"
export SPECIES_METADATA="${METADATA_DIR}/species_list_beta2_signature.tsv"

export OUT_DIR="${PROJECT_DIR}/results/04_evolutionary_analysis"
export INPUT_DIR="${OUT_DIR}/input"
export ORTHOFINDER_DIR="${INPUT_DIR}/orthofinder"
export KAKS_DIR="${INPUT_DIR}/kaks"
export BETA2_DIR="${INPUT_DIR}/beta2"
export PHYLOGENY_DIR="${INPUT_DIR}/phylogeny"
export LOG_DIR="${OUT_DIR}/logs"

export CROSS_SPECIES_DBD_FASTA="${BETA2_DIR}/cross_species_erf_ap2_domains_aligned.fa"
export CROSS_SPECIES_BETA2="${BETA2_DIR}/cross_species_beta2_signature.tsv"

mkdir -p "${OUT_DIR}" "${INPUT_DIR}" "${ORTHOFINDER_DIR}" "${KAKS_DIR}" "${BETA2_DIR}" "${PHYLOGENY_DIR}" "${LOG_DIR}"

echo "PROJECT_DIR=${PROJECT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
