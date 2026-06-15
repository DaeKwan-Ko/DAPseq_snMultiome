#!/usr/bin/env bash
# Convert JGI/MACS-style peak files to BED files for ChIPseeker.
# Duplicated intervals are collapsed by keeping the row with the highest score.
# Usage: bash 04_prepare_peak_beds_for_chipseeker.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

mkdir -p "${OUT_DIR}/peaks_bed_for_chipseeker"

if [[ ! -d "${PEAK_DIR}" ]]; then
    echo "ERROR: PEAK_DIR does not exist: ${PEAK_DIR}" >&2
    exit 1
fi

# Accept common peak file suffixes. Adjust if JGI used a different extension.
shopt -s nullglob
peak_files=("${PEAK_DIR}"/*.xls "${PEAK_DIR}"/*.bed "${PEAK_DIR}"/*.narrowPeak)

if [[ ${#peak_files[@]} -eq 0 ]]; then
    echo "ERROR: no peak files found in ${PEAK_DIR}" >&2
    exit 1
fi

for peak_file in "${peak_files[@]}"; do
    base="$(basename "${peak_file}")"
    gene="$(echo "${base}" | sed 's/[_].*//' | sed 's/\.xls$//' | sed 's/\.bed$//' | sed 's/\.narrowPeak$//')"
    out="${OUT_DIR}/peaks_bed_for_chipseeker/${gene}.bed"

    # Output BED5: chrom, start, end, name, score.
    # Skip comment/header lines. Keep only rows with numeric start/end.
    # If no score-like column exists, use 1.
    awk -v OFS='\t' -v gene="${gene}" '
        BEGIN {n=0}
        /^#/ {next}
        NR==1 && ($1 ~ /chr|chrom|name/i || $2 ~ /start/i) {next}
        NF >= 3 && $2 ~ /^[0-9]+$/ && $3 ~ /^[0-9]+$/ {
            score=1
            if (NF >= 5 && $5 ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) score=$5
            else if (NF >= 7 && $7 ~ /^[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$/) score=$7
            name=gene"_peak"(++n)
            print $1,$2,$3,name,score
        }
    ' "${peak_file}" \
    | sort -k1,1 -k2,2n -k3,3n -k5,5nr \
    | awk -v OFS='\t' '!seen[$1 FS $2 FS $3]++ {print}' \
    > "${out}"

    echo "Wrote ${out}"
done

echo "BED files for ChIPseeker are in: ${OUT_DIR}/peaks_bed_for_chipseeker"
