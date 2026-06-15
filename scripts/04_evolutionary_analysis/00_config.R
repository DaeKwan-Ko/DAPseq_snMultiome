# 00_config.R
# Central R configuration for evolutionary analyses.

suppressPackageStartupMessages({
  library(fs)
})

PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)

METADATA_DIR <- file.path(PROJECT_DIR, "metadata")
ERF_CLADE_METADATA <- file.path(METADATA_DIR, "erf_clade_annotation.tsv")
SPECIES_METADATA <- file.path(METADATA_DIR, "species_list_beta2_signature.tsv")

# Upstream outputs from 03_erf_phylogeny_motifs
ERF_META_PREPARED <- file.path(PROJECT_DIR, "results", "03_erf_phylogeny_motifs", "erf_metadata_prepared.tsv")
SORGHUM_BETA2_SIGNATURE <- file.path(PROJECT_DIR, "results", "03_erf_phylogeny_motifs", "erf_beta2_signature.tsv")

# Evolutionary analysis input directories
OUT_DIR <- file.path(PROJECT_DIR, "results", "04_evolutionary_analysis")
INPUT_DIR <- file.path(OUT_DIR, "input")
ORTHOFINDER_DIR <- file.path(INPUT_DIR, "orthofinder")
KAKS_DIR <- file.path(INPUT_DIR, "kaks")
BETA2_DIR <- file.path(INPUT_DIR, "beta2")
PHYLOGENY_DIR <- file.path(INPUT_DIR, "phylogeny")

# Expected input files
KAKS_RESULTS <- file.path(KAKS_DIR, "ap2_dbd_kaks_results.tsv")
CROSS_SPECIES_BETA2 <- file.path(BETA2_DIR, "cross_species_beta2_signature.tsv")
CROSS_SPECIES_DBD_FASTA <- file.path(BETA2_DIR, "cross_species_erf_ap2_domains_aligned.fa")

# Output directories
dir_create(OUT_DIR, recurse = TRUE)
dir_create(INPUT_DIR, recurse = TRUE)
dir_create(ORTHOFINDER_DIR, recurse = TRUE)
dir_create(KAKS_DIR, recurse = TRUE)
dir_create(BETA2_DIR, recurse = TRUE)
dir_create(PHYLOGENY_DIR, recurse = TRUE)

# Labels used throughout the manuscript
CLADE_LABELS <- c(A = "ECA", B = "ECB")

message("PROJECT_DIR: ", PROJECT_DIR)
message("OUT_DIR: ", OUT_DIR)
