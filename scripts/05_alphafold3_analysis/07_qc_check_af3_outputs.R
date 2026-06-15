#!/usr/bin/env Rscript

# 07_qc_check_af3_outputs.R
# Quick QC for AlphaFold3 analysis outputs.

source("scripts/05_alphafold3_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

check_file <- function(path) {
  data.frame(
    file = basename(path),
    path = path,
    exists = file.exists(path),
    size_bytes = ifelse(file.exists(path), file.info(path)$size, NA_real_)
  )
}

paths <- c(
  AF3_DESIGNS,
  MANIFEST,
  CONFIDENCE_SUMMARY,
  file.path(OUT_DIR, "af3_atom_table_manifest.tsv"),
  DNA_RMSD,
  CONTACTS,
  BETA2_DISTANCES,
  file.path(OUT_DIR, "af3_dna_rmsd_summary.tsv"),
  file.path(OUT_DIR, "af3_base_contact_summary.tsv"),
  file.path(OUT_DIR, "af3_beta2_distance_summary.tsv")
)

qc <- bind_rows(lapply(paths, check_file))
write_tsv(qc, file.path(OUT_DIR, "af3_analysis_output_qc.tsv"))
print(qc)

if (file.exists(DNA_RMSD)) {
  rmsd <- read_tsv(DNA_RMSD, show_col_types = FALSE)
  message("DNA RMSD rows: ", nrow(rmsd))
  if ("design_id" %in% colnames(rmsd)) print(table(rmsd$design_id))
}

if (file.exists(CONTACTS)) {
  contacts <- read_tsv(CONTACTS, show_col_types = FALSE)
  message("Contact rows: ", nrow(contacts))
}

if (file.exists(BETA2_DISTANCES)) {
  beta2 <- read_tsv(BETA2_DISTANCES, show_col_types = FALSE)
  message("Beta2 distance rows: ", nrow(beta2))
}

message("Saved QC table: ", file.path(OUT_DIR, "af3_analysis_output_qc.tsv"))
