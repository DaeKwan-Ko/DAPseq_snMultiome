#!/usr/bin/env bash
# Configuration file for DAP-seq post-processing scripts.
# Edit all paths below before running the workflow.

set -euo pipefail

# Project root for local or HPCC analysis.
PROJECT_DIR="/mnt/scratch/$USER/glbrc/dapseq"

# JGI-provided DAP-seq processed data directory.
JGI_DIR="${PROJECT_DIR}/jgi_provided/DAP_analysis"

# Directory containing JGI-provided BAM files and BAM index files.
ALIGN_DIR="${JGI_DIR}/Alignments"

# Directory containing JGI-provided peak files.
# Change this if your peak files are in a different subdirectory.
PEAK_DIR="${JGI_DIR}/Peaks"

# Output directory for this workflow.
OUT_DIR="${PROJECT_DIR}/analysis/01_dapseq_processing"

# Metadata table with at least two columns: geneID and tf_family.
# Example: metadata/tf_list_142.tsv
TF_METADATA="${PROJECT_DIR}/metadata/tf_list_142.tsv"

# Sorghum reference genome and annotation.
GENOME_FASTA="${PROJECT_DIR}/reference/Sbicolor_454_v3.0.fa"
GENOME_GFF3="${PROJECT_DIR}/reference/Sbicolor_454_v3.1.1.gene.gff3"

# Effective genome size used for MACS/JGI reporting if needed.
# This is not used by all scripts.
GENOME_SIZE="732200000"

# deepTools settings.
BIN_SIZE=20
THREADS=16

# Names of the major TF families used for family-specific correlation analysis.
FAMILIES="ERF MYB WRKY"

mkdir -p "${OUT_DIR}"
