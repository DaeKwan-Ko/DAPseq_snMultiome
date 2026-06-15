#!/usr/bin/env Rscript

# 10_qc_check_snmultiome_outputs.R
# Quick QC checks for snMultiome processing outputs.

source("scripts/06_snmultiome_processing/00_config.R")

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
  SAMPLE_METADATA,
  file.path(OUT_DIR, "snmultiome_sample_objects_raw.rds"),
  MERGED_OBJECT,
  WNN_OBJECT,
  ANNOTATED_OBJECT,
  GENE_ACTIVITY_OBJECT,
  file.path(OUT_DIR, "joint_qc_summary_by_sample.tsv"),
  file.path(OUT_DIR, "celltype_annotation_metadata.tsv"),
  file.path(OUT_DIR, "celltype_condition_nucleus_counts.tsv"),
  file.path(OUT_DIR, "rna_by_celltype_condition.tsv"),
  file.path(OUT_DIR, "activity_by_celltype_condition.tsv"),
  file.path(OUT_DIR, "erf_rna_activity_summary_by_celltype_condition.tsv"),
  file.path(OUT_DIR, "erf_rna_activity_log2fc_by_celltype.tsv"),
  file.path(OUT_DIR, "target_gene_rna_activity_summary_by_celltype_condition.tsv")
)

qc <- bind_rows(lapply(paths, check_file))
write_tsv(qc, file.path(OUT_DIR, "snmultiome_processing_qc.tsv"))
print(qc)

counts_file <- file.path(OUT_DIR, "celltype_condition_nucleus_counts.tsv")
if (file.exists(counts_file)) {
  counts <- read_tsv(counts_file, show_col_types = FALSE)
  message("Celltype × condition nuclei counts:")
  print(counts)
}

joint_qc_file <- file.path(OUT_DIR, "joint_qc_summary_by_sample.tsv")
if (file.exists(joint_qc_file)) {
  joint_qc <- read_tsv(joint_qc_file, show_col_types = FALSE)
  message("Joint QC summary:")
  print(joint_qc)
}

message("Saved QC table: ", file.path(OUT_DIR, "snmultiome_processing_qc.tsv"))
