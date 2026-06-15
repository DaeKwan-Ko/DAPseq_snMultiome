#!/usr/bin/env Rscript

# 07_qc_check_annotation_outputs.R
# Quick QC checks for DAP-seq annotation outputs.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

files_to_check <- c(
  "all_peak_annotations.tsv",
  "peak_annotation_summary_by_tf.tsv",
  "dapseq_candidate_targets_long.tsv",
  "dapseq_candidate_targets_by_tf.tsv",
  "dapseq_target_count_by_tf.tsv",
  "erf_target_matrix.tsv",
  "erf_target_clade_counts.tsv",
  "erf_clade_biased_recurrent_targets.tsv"
)

qc <- lapply(files_to_check, function(f) {
  path <- file.path(OUT_DIR, f)
  if (!file.exists(path)) {
    return(data.frame(file = f, exists = FALSE, n_rows = NA_integer_, n_cols = NA_integer_))
  }
  x <- read_tsv(path, show_col_types = FALSE)
  data.frame(file = f, exists = TRUE, n_rows = nrow(x), n_cols = ncol(x))
}) %>%
  bind_rows()

write_tsv(qc, file.path(OUT_DIR, "annotation_output_qc.tsv"))

print(qc)

message("Saved QC table: ", file.path(OUT_DIR, "annotation_output_qc.tsv"))
