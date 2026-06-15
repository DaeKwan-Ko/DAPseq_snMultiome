#!/usr/bin/env bash
# Generate BAM file lists and labels for deepTools multiBamSummary.
# Run from any directory after editing 00_config.sh.
# Usage: bash 01_generate_bam_filelists.sh

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${SCRIPT_DIR}/00_config.sh"

mkdir -p "${OUT_DIR}/filelists"

if [[ ! -d "${ALIGN_DIR}" ]]; then
    echo "ERROR: ALIGN_DIR does not exist: ${ALIGN_DIR}" >&2
    exit 1
fi

if [[ ! -f "${TF_METADATA}" ]]; then
    echo "ERROR: TF metadata file does not exist: ${TF_METADATA}" >&2
    exit 1
fi

# Metadata should have columns named geneID and tf_family.
# If your file uses label/tf, rename them before running or edit the awk commands below.
awk 'NR==1 {for (i=1; i<=NF; i++) {if ($i=="geneID") gid=i; if ($i=="tf_family") fam=i} next} {print $gid"\t"$fam}' \
    "${TF_METADATA}" > "${OUT_DIR}/filelists/tf_family_map.tsv"

# All BAMs sorted by sorghum gene ID component.
find "${ALIGN_DIR}" -maxdepth 1 -type f -name 'Sobic*.bam' | sort -t. -k2 > "${OUT_DIR}/filelists/all_bams.txt"

# Labels are gene IDs inferred from BAM filenames before the first underscore.
sed 's#.*/##' "${OUT_DIR}/filelists/all_bams.txt" | sed 's/[_].*//' > "${OUT_DIR}/filelists/all_labels.txt"

# Space-delimited versions for easy use in shell commands.
paste -sd ' ' "${OUT_DIR}/filelists/all_bams.txt" > "${OUT_DIR}/filelists/all_bams.space.txt"
paste -sd ' ' "${OUT_DIR}/filelists/all_labels.txt" > "${OUT_DIR}/filelists/all_labels.space.txt"

# Family-specific file lists.
for family in ${FAMILIES}; do
    family_lc="$(echo "${family}" | tr '[:upper:]' '[:lower:]')"
    family_genes="${OUT_DIR}/filelists/${family_lc}_genes.txt"
    awk -v f="${family}" '$2==f {print $1}' "${OUT_DIR}/filelists/tf_family_map.tsv" | sort -t. -k2 > "${family_genes}"

    : > "${OUT_DIR}/filelists/${family_lc}_bams.txt"
    : > "${OUT_DIR}/filelists/${family_lc}_labels.txt"

    while read -r gene; do
        bam="$(find "${ALIGN_DIR}" -maxdepth 1 -type f -name "${gene}_*.bam" | head -n 1 || true)"
        if [[ -n "${bam}" ]]; then
            echo "${bam}" >> "${OUT_DIR}/filelists/${family_lc}_bams.txt"
            echo "${gene}" >> "${OUT_DIR}/filelists/${family_lc}_labels.txt"
        else
            echo "WARNING: no BAM found for ${gene}" >&2
        fi
    done < "${family_genes}"

    paste -sd ' ' "${OUT_DIR}/filelists/${family_lc}_bams.txt" > "${OUT_DIR}/filelists/${family_lc}_bams.space.txt"
    paste -sd ' ' "${OUT_DIR}/filelists/${family_lc}_labels.txt" > "${OUT_DIR}/filelists/${family_lc}_labels.space.txt"
done

echo "Generated BAM file lists in: ${OUT_DIR}/filelists"
