#!/usr/bin/env bash
set -euo pipefail

# 07_make_ceqlogo_commands.sh
# Generate ceqlogo commands for MEME motif files.
#
# This script creates a command file:
#   results/03_erf_phylogeny_motifs/ceqlogo_commands.sh
#
# Review the command file before running it. Add reverse-complement commands
# manually for motifs that need reverse-complement display.

source scripts/03_erf_phylogeny_motifs/00_config.sh

COMMAND_FILE="${OUT_DIR}/ceqlogo_commands.sh"
: > "${COMMAND_FILE}"

if [[ ! -d "${MEME_DIR}" ]]; then
  echo "MEME_DIR not found: ${MEME_DIR}" >&2
  exit 0
fi

find "${MEME_DIR}" -type f \( -name "meme.txt" -o -name "*.meme" \) | sort | while read -r meme_file; do
  gene_id=$(echo "${meme_file}" | grep -oE 'Sobic\.[0-9]{3}G[0-9]+' | head -n 1 || true)
  if [[ -z "${gene_id}" ]]; then
    gene_id=$(basename "$(dirname "${meme_file}")")
  fi

  out_png="${LOGO_DIR}/${gene_id}_motif1_logo.png"
  echo "${CEQLOGO_BIN} -i \"${meme_file}\" -m 1 -o \"${out_png}\" -f PNG" >> "${COMMAND_FILE}"
done

chmod +x "${COMMAND_FILE}"

echo "Saved ceqlogo command file: ${COMMAND_FILE}"
echo "Review it before running. To reverse-complement a motif, add -r to that command."
