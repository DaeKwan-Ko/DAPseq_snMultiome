#!/usr/bin/env bash

# 00_config.sh
# Shell configuration for Cell Ranger ARC helper scripts.

export PROJECT_DIR="${PROJECT_DIR:-$(pwd)}"

export METADATA_DIR="${PROJECT_DIR}/metadata"
export OUT_DIR="${PROJECT_DIR}/results/06_snmultiome_processing"
export INPUT_DIR="${OUT_DIR}/input"
export FASTQ_DIR="${INPUT_DIR}/fastq"
export CELLRANGER_ARC_DIR="${INPUT_DIR}/cellranger_arc"

# Reference-building inputs
export GENOME_FASTA="${PROJECT_DIR}/reference/Sbicolor_454_v3.0.1.fa"
export GENOME_GFF3="${PROJECT_DIR}/reference/Sbicolor_454_v3.1.1.gene.gff3"
export ARC_REFERENCE_DIR="${PROJECT_DIR}/reference/cellranger_arc_sorghum_v3.1.1"

# Cell Ranger ARC executable
export CELLRANGER_ARC_BIN="${CELLRANGER_ARC_BIN:-cellranger-arc}"

# Sample metadata
export SAMPLE_METADATA="${METADATA_DIR}/sample_metadata_snmultiome.tsv"

mkdir -p "${OUT_DIR}" "${INPUT_DIR}" "${FASTQ_DIR}" "${CELLRANGER_ARC_DIR}"

echo "PROJECT_DIR=${PROJECT_DIR}"
echo "OUT_DIR=${OUT_DIR}"
