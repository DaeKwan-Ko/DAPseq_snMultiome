#!/usr/bin/env Rscript

# 01_prepare_species_metadata.R
# Prepare species-level metadata for cross-lineage ERF β2 signature analyses.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

if (!file.exists(SPECIES_METADATA)) {
  stop("Species metadata file not found: ", SPECIES_METADATA)
}

species_meta <- read_tsv(SPECIES_METADATA, show_col_types = FALSE)

required_cols <- c("lineage_order", "species", "taxonomy", "major_group", "total_erfs_analyzed")
missing_cols <- setdiff(required_cols, colnames(species_meta))
if (length(missing_cols) > 0) {
  stop("Missing required columns in species metadata: ", paste(missing_cols, collapse = ", "))
}

species_meta <- species_meta %>%
  mutate(
    species_label = str_replace_all(species, " ", "_"),
    species_display = species,
    lineage_order = as.integer(lineage_order),
    total_erfs_analyzed = as.integer(total_erfs_analyzed)
  ) %>%
  arrange(lineage_order)

out_file <- file.path(OUT_DIR, "species_metadata_prepared.tsv")
write_tsv(species_meta, out_file)

message("Saved species metadata: ", out_file)
message("Number of species: ", nrow(species_meta))
message("Total ERFs expected across species: ", sum(species_meta$total_erfs_analyzed, na.rm = TRUE))
message("Done.")
