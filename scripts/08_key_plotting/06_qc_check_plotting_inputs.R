#!/usr/bin/env Rscript

# 06_qc_check_plotting_inputs.R
# Check whether key plotting inputs and outputs exist.

source("scripts/08_key_plotting/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

check_file <- function(path, category) {
  data.frame(
    category = category,
    file = basename(path),
    path = path,
    exists = file.exists(path),
    size_bytes = ifelse(file.exists(path), file.info(path)$size, NA_real_)
  )
}

inputs <- c(
  PEAK_ANNOT_SUMMARY,
  ERF_TARGET_CLADE_COUNTS,
  ERF_CLADE_BIASED_TARGETS,
  ERF_PHYLO_MOTIF_TABLE,
  ERF_TREE_NWK,
  MEME_TOP_MOTIFS,
  KAKS_PAIRS,
  KAKS_SUMMARY,
  KAKS_TESTS,
  BETA2_LINEAGE_SUMMARY,
  CELLTYPE_COUNTS,
  ERF_RNA_ACTIVITY_FC,
  TARGET_RNA_ACTIVITY_FC,
  NETWORK_EDGES,
  DRIVER_RANKING,
  WEIGHTED_OUT_DEGREE
)

plots <- list.files(PLOT_DIR, pattern = "\\.pdf$", full.names = TRUE)
source_data <- list.files(SOURCE_DATA_DIR, pattern = "\\.tsv$", full.names = TRUE)

qc <- bind_rows(
  bind_rows(lapply(inputs, check_file, category = "input")),
  bind_rows(lapply(plots, check_file, category = "plot")),
  bind_rows(lapply(source_data, check_file, category = "source_data"))
)

write_tsv(qc, file.path(OUT_DIR, "plotting_input_qc.tsv"))
print(qc)

message("Saved plotting QC table: ", file.path(OUT_DIR, "plotting_input_qc.tsv"))
