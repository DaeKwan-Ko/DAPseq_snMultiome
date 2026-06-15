#!/usr/bin/env Rscript

# 05_prepare_source_data_tables.R
# Prepare compact source-data tables for key figure-style quantitative panels.

source("scripts/08_key_plotting/00_config.R")
source("scripts/08_key_plotting/plotting_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

copy_if_exists <- function(input, output_name) {
  if (file.exists(input)) {
    x <- read_tsv(input, show_col_types = FALSE)
    out <- file.path(SOURCE_DATA_DIR, output_name)
    write_tsv(x, out)
    message("Wrote source data: ", out)
    return(tibble(source_file = input, output_file = out, n_rows = nrow(x), n_cols = ncol(x)))
  }
  tibble(source_file = input, output_file = NA_character_, n_rows = NA_integer_, n_cols = NA_integer_)
}

manifest <- bind_rows(
  copy_if_exists(PEAK_ANNOT_SUMMARY, "source_peak_annotation_summary_by_tf.tsv"),
  copy_if_exists(ERF_TARGET_CLADE_COUNTS, "source_erf_target_clade_counts.tsv"),
  copy_if_exists(ERF_PHYLO_MOTIF_TABLE, "source_erf_phylogeny_motif_annotation_table.tsv"),
  copy_if_exists(KAKS_PAIRS, "source_kaks_ortholog_pairs_annotated.tsv"),
  copy_if_exists(BETA2_LINEAGE_SUMMARY, "source_cross_lineage_beta2_summary.tsv"),
  copy_if_exists(CELLTYPE_COUNTS, "source_celltype_condition_nucleus_counts.tsv"),
  copy_if_exists(ERF_RNA_ACTIVITY_FC, "source_erf_rna_activity_log2fc_by_celltype.tsv"),
  copy_if_exists(TARGET_RNA_ACTIVITY_FC, "source_target_rna_activity_log2fc_by_celltype.tsv"),
  copy_if_exists(NETWORK_EDGES, "source_cortex1_erf_target_edges_scored.tsv"),
  copy_if_exists(DRIVER_RANKING, "source_cortex1_erf_driver_ranking.tsv")
)

write_tsv(manifest, file.path(SOURCE_DATA_DIR, "source_data_manifest.tsv"))

message("Saved source-data manifest: ", file.path(SOURCE_DATA_DIR, "source_data_manifest.tsv"))
