#!/usr/bin/env bash
set -euo pipefail

# 02_run_cellranger_arc_count.sh
# Run Cell Ranger ARC count for each sample in sample_metadata_snmultiome.tsv.
#
# Run from repository root:
#   bash scripts/06_snmultiome_processing/02_run_cellranger_arc_count.sh
#
# This script assumes FASTQ files are already downloaded and organized under:
#   results/06_snmultiome_processing/input/fastq/<sample_id>/
#
# Large FASTQ and Cell Ranger ARC outputs should not be committed to GitHub.

source scripts/06_snmultiome_processing/00_config.sh

if [[ ! -s "${SAMPLE_METADATA}" ]]; then
  echo "ERROR: sample metadata not found: ${SAMPLE_METADATA}" >&2
  exit 1
fi

if [[ ! -d "${ARC_REFERENCE_DIR}" ]]; then
  echo "ERROR: Cell Ranger ARC reference not found: ${ARC_REFERENCE_DIR}" >&2
  exit 1
fi

tail -n +2 "${SAMPLE_METADATA}" | while IFS=$'\t' read -r sample_id condition treatment rest; do
  if [[ -z "${sample_id}" ]]; then
    continue
  fi

  SAMPLE_FASTQ_DIR="${FASTQ_DIR}/${sample_id}"
  SAMPLE_OUT_DIR="${CELLRANGER_ARC_DIR}/${sample_id}"

  if [[ ! -d "${SAMPLE_FASTQ_DIR}" ]]; then
    echo "WARNING: FASTQ directory not found for ${sample_id}: ${SAMPLE_FASTQ_DIR}; skipping." >&2
    continue
  fi

  mkdir -p "${CELLRANGER_ARC_DIR}"
  cd "${CELLRANGER_ARC_DIR}"

  echo "Running Cell Ranger ARC count for ${sample_id}"

  "${CELLRANGER_ARC_BIN}" count \
    --id="${sample_id}" \
    --reference="${ARC_REFERENCE_DIR}" \
    --libraries="${SAMPLE_FASTQ_DIR}/libraries.csv" \
    --localcores=16 \
    --localmem=120

  echo "Finished ${sample_id}: ${SAMPLE_OUT_DIR}"
  cd "${PROJECT_DIR}"
done
