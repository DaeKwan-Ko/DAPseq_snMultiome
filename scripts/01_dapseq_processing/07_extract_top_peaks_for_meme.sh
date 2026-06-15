#!/usr/bin/env bash
# Extract top-ranked peak sequences for MEME motif discovery.
# The script writes two FASTA sets per TF: full peak regions and summit-centered regions.
# Usage: bash 07_extract_top_peaks_for_meme.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

BED_DIR="${OUT_DIR}/peaks_bed_for_chipseeker"
MEME_IN_DIR="${OUT_DIR}/meme_inputs"
TOP_N=100
SUMMIT_FLANK=30

mkdir -p "${MEME_IN_DIR}/full_peaks" "${MEME_IN_DIR}/summit_pm${SUMMIT_FLANK}bp"

if [[ ! -f "${GENOME_FASTA}" ]]; then
    echo "ERROR: GENOME_FASTA not found: ${GENOME_FASTA}" >&2
    exit 1
fi

if [[ ! -d "${BED_DIR}" ]]; then
    echo "ERROR: BED_DIR does not exist. Run 04_prepare_peak_beds_for_chipseeker.sh first." >&2
    exit 1
fi

# Requires bedtools getfasta.
# module load BEDTools

for bed in "${BED_DIR}"/*.bed; do
    gene="$(basename "${bed}" .bed)"
    top_bed="${MEME_IN_DIR}/${gene}.top${TOP_N}.bed"
    summit_bed="${MEME_IN_DIR}/${gene}.top${TOP_N}.summit_pm${SUMMIT_FLANK}bp.bed"

    # Sort by score descending and retain top N peaks.
    sort -k5,5nr "${bed}" | head -n "${TOP_N}" > "${top_bed}"

    # Approximate summit as midpoint if summit coordinates are not available in the BED5 file.
    awk -v OFS='\t' -v flank="${SUMMIT_FLANK}" '{mid=int(($2+$3)/2); start=mid-flank; if (start<0) start=0; end=mid+flank; print $1,start,end,$4,$5}' "${top_bed}" > "${summit_bed}"

    bedtools getfasta -fi "${GENOME_FASTA}" -bed "${top_bed}" -name -fo "${MEME_IN_DIR}/full_peaks/${gene}.top${TOP_N}.fa"
    bedtools getfasta -fi "${GENOME_FASTA}" -bed "${summit_bed}" -name -fo "${MEME_IN_DIR}/summit_pm${SUMMIT_FLANK}bp/${gene}.top${TOP_N}.summit_pm${SUMMIT_FLANK}bp.fa"
done

echo "MEME input FASTA files written to: ${MEME_IN_DIR}"
