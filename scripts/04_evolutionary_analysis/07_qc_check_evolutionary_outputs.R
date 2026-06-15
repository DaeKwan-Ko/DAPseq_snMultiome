#!/usr/bin/env Rscript

# 07_qc_check_evolutionary_outputs.R
# Quick QC check for evolutionary analysis outputs.

source("scripts/04_evolutionary_analysis/00_config.R")

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
  file.path(OUT_DIR, "species_metadata_prepared.tsv"),
  file.path(OUT_DIR, "ortholog_pairs_long.tsv"),
  file.path(OUT_DIR, "kaks_ortholog_pairs_annotated.tsv"),
  file.path(OUT_DIR, "kaks_summary_by_clade.tsv"),
  file.path(OUT_DIR, "kaks_wilcoxon_tests.tsv"),
  file.path(OUT_DIR, "cross_species_beta2_signature_long.tsv"),
  file.path(OUT_DIR, "cross_species_beta2_signature_counts.tsv"),
  file.path(OUT_DIR, "cross_lineage_beta2_summary_for_plotting.tsv"),
  file.path(OUT_DIR, "combined_erf_phylogeny_metadata.tsv")
)

qc <- bind_rows(lapply(paths, check_file))
write_tsv(qc, file.path(OUT_DIR, "evolutionary_analysis_output_qc.tsv"))
print(qc)

# Extra helpful checks
kaks_file <- file.path(OUT_DIR, "kaks_ortholog_pairs_annotated.tsv")
if (file.exists(kaks_file)) {
  kaks <- read_tsv(kaks_file, show_col_types = FALSE)
  message("Usable Ka/Ks pairs: ", nrow(kaks))
  if (all(c("comparison", "clade") %in% colnames(kaks))) {
    print(table(kaks$comparison, kaks$clade, useNA = "ifany"))
  }
}

beta2_file <- file.path(OUT_DIR, "cross_lineage_beta2_summary_for_plotting.tsv")
if (file.exists(beta2_file)) {
  beta2 <- read_tsv(beta2_file, show_col_types = FALSE)
  message("Cross-lineage beta2 species: ", nrow(beta2))
  if ("n_total_classified" %in% colnames(beta2)) {
    message("Total ERFs classified: ", sum(beta2$n_total_classified, na.rm = TRUE))
  }
}

message("Saved QC table: ", file.path(OUT_DIR, "evolutionary_analysis_output_qc.tsv"))
