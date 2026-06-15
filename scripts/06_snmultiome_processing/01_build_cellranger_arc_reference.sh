#!/usr/bin/env bash
set -euo pipefail

# 01_build_cellranger_arc_reference.sh
# Build a Cell Ranger ARC reference for Sorghum bicolor.
#
# Run from repository root:
#   bash scripts/06_snmultiome_processing/01_build_cellranger_arc_reference.sh
#
# Edit scripts/06_snmultiome_processing/00_config.sh before running.
#
# Note:
# Cell Ranger ARC expects a genome FASTA and GTF. If your annotation is GFF3,
# convert it to a Cell Ranger-compatible GTF before running this script.

source scripts/06_snmultiome_processing/00_config.sh

GTF_FILE="${PROJECT_DIR}/reference/Sbicolor_454_v3.1.1.gene.gtf"

if [[ ! -s "${GENOME_FASTA}" ]]; then
  echo "ERROR: GENOME_FASTA not found: ${GENOME_FASTA}" >&2
  exit 1
fi

if [[ ! -s "${GTF_FILE}" ]]; then
  echo "ERROR: GTF file not found: ${GTF_FILE}" >&2
  echo "Convert the GFF3 annotation to GTF first, then rerun." >&2
  exit 1
fi

mkdir -p "$(dirname "${ARC_REFERENCE_DIR}")"

"${CELLRANGER_ARC_BIN}" mkref \
  --config <(cat <<EOF
{
  "genome": "Sorghum_bicolor_v3_1_1",
  "fasta": "${GENOME_FASTA}",
  "genes": "${GTF_FILE}"
}
EOF
) \
  --nthreads=16 \
  --memgb=120 \
  --ref-version="Sorghum_bicolor_v3.1.1"

echo "Cell Ranger ARC reference created at: ${ARC_REFERENCE_DIR}"
