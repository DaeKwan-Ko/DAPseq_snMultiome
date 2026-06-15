#!/usr/bin/env Rscript

# 06_prepare_combined_erf_phylogeny_metadata.R
# Prepare metadata table for a combined cross-lineage AP2/ERF DBD phylogeny.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

species_file <- file.path(OUT_DIR, "species_metadata_prepared.tsv")
beta2_long_file <- file.path(OUT_DIR, "cross_species_beta2_signature_long.tsv")

if (!file.exists(species_file)) {
  stop("Prepared species metadata not found. Run 01_prepare_species_metadata.R first.")
}

if (!file.exists(beta2_long_file)) {
  message("Cross-species beta2 long table not found. Run 05_build_cross_species_beta2_summary.R first.")
  quit(save = "no", status = 0)
}

species_meta <- read_tsv(species_file, show_col_types = FALSE)
beta2 <- read_tsv(beta2_long_file, show_col_types = FALSE)

phylo_meta <- beta2 %>%
  mutate(species = str_replace_all(species, "_", " ")) %>%
  left_join(species_meta, by = "species") %>%
  mutate(
    tip_label = paste(str_replace_all(species, " ", "_"), geneID, sep = "|"),
    beta2_plot_class = factor(beta2_signature_class, levels = c("AA", "WV", "neither")),
    taxonomy_plot_order = lineage_order
  ) %>%
  arrange(lineage_order, species, beta2_signature_class, geneID)

out_file <- file.path(OUT_DIR, "combined_erf_phylogeny_metadata.tsv")
write_tsv(phylo_meta, out_file)

message("Saved combined phylogeny metadata: ", out_file)
message("Rows: ", nrow(phylo_meta))
message("Done.")
