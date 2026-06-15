# 00_config.R
# Central R configuration for AlphaFold3 analysis summary scripts.

suppressPackageStartupMessages({
  library(fs)
})

PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)

METADATA_DIR <- file.path(PROJECT_DIR, "metadata")
OUT_DIR <- file.path(PROJECT_DIR, "results", "05_alphafold3_analysis")
INPUT_DIR <- file.path(OUT_DIR, "input")
AF3_OUTPUT_DIR <- file.path(INPUT_DIR, "af3_outputs")

AF3_DESIGNS <- file.path(METADATA_DIR, "af3_designs.tsv")
MOTIF_POSITION_MAP <- file.path(METADATA_DIR, "af3_motif_position_map.tsv")
BETA2_RESIDUE_MAP <- file.path(METADATA_DIR, "af3_beta2_residue_map.tsv")

MANIFEST <- file.path(OUT_DIR, "af3_output_manifest.tsv")
CONFIDENCE_SUMMARY <- file.path(OUT_DIR, "af3_confidence_summary.tsv")
DNA_RMSD <- file.path(OUT_DIR, "dna_rmsd_by_model.tsv")
CONTACTS <- file.path(OUT_DIR, "protein_dna_base_contacts.tsv")
BETA2_DISTANCES <- file.path(OUT_DIR, "beta2_to_dna_base_min_distances.tsv")

dir_create(OUT_DIR, recurse = TRUE)
dir_create(file.path(OUT_DIR, "plots"), recurse = TRUE)

message("PROJECT_DIR: ", PROJECT_DIR)
message("OUT_DIR: ", OUT_DIR)
