#!/usr/bin/env Rscript

# 01_prepare_erf_metadata.R
# Prepare a clean metadata table for the 49 DAP-seq-profiled sorghum ERFs.

source("scripts/03_erf_phylogeny_motifs/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

if (!file.exists(TF_METADATA)) {
  stop("TF metadata file not found: ", TF_METADATA)
}
if (!file.exists(ERF_CLADE_METADATA)) {
  stop("ERF clade metadata file not found: ", ERF_CLADE_METADATA)
}

tf <- read_tsv(TF_METADATA, show_col_types = FALSE)

# Accept older headers if needed.
if (all(c("label", "tf") %in% colnames(tf)) && !("geneID" %in% colnames(tf))) {
  tf <- tf %>% rename(geneID = label, tf_family = tf)
}

clade <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE)

if (!all(c("geneID", "clade") %in% colnames(clade))) {
  stop("ERF clade metadata must contain columns: geneID, clade")
}

salinity <- NULL
if (file.exists(SALINITY_ERF_METADATA)) {
  salinity <- read_tsv(SALINITY_ERF_METADATA, show_col_types = FALSE)
  if (!all(c("geneID", "clade") %in% colnames(salinity))) {
    warning("Salinity ERF metadata does not contain expected columns geneID and clade; ignoring.")
    salinity <- NULL
  }
}

erf_meta <- tf %>%
  filter(tf_family == "ERF") %>%
  select(geneID, tf_family) %>%
  left_join(clade %>% select(geneID, clade), by = "geneID") %>%
  mutate(
    clade_full = case_when(
      clade == "A" ~ "ECA",
      clade == "B" ~ "ECB",
      TRUE ~ NA_character_
    ),
    salinity_linked = if (!is.null(salinity)) geneID %in% salinity$geneID else NA
  ) %>%
  arrange(clade, geneID)

if (any(is.na(erf_meta$clade))) {
  warning("Some ERFs are missing clade assignments.")
}

out_file <- file.path(OUT_DIR, "erf_metadata_prepared.tsv")
write_tsv(erf_meta, out_file)

message("Saved ERF metadata: ", out_file)
message("Number of ERFs: ", nrow(erf_meta))
print(table(erf_meta$clade, useNA = "ifany"))
message("Done.")
