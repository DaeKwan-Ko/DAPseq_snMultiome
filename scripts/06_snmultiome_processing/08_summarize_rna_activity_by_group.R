#!/usr/bin/env Rscript

# 08_summarize_rna_activity_by_group.R
# Summarize RNA expression and ATAC-derived gene activity by celltype × condition.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(Matrix)
})

if (!file.exists(GENE_ACTIVITY_OBJECT)) {
  stop("Object with gene activity assay not found. Run 07_calculate_gene_activity.R first: ", GENE_ACTIVITY_OBJECT)
}

obj <- readRDS(GENE_ACTIVITY_OBJECT)

if (!all(c("celltype", "condition") %in% colnames(obj@meta.data))) {
  stop("Object metadata must contain celltype and condition columns.")
}

obj$celltype_condition <- paste(obj$celltype, obj$condition, sep = "__")

summarize_assay <- function(obj, assay_name, slot_name = "data") {
  DefaultAssay(obj) <- assay_name
  mat <- GetAssayData(obj, assay = assay_name, slot = slot_name)

  groups <- obj$celltype_condition
  group_levels <- sort(unique(groups))

  out <- lapply(group_levels, function(g) {
    cells <- colnames(obj)[groups == g]
    if (length(cells) == 0) return(NULL)

    sub <- mat[, cells, drop = FALSE]
    avg <- Matrix::rowMeans(sub)
    pct <- Matrix::rowMeans(sub > 0) * 100

    parts <- strsplit(g, "__", fixed = TRUE)[[1]]
    data.frame(
      geneID = rownames(mat),
      celltype = parts[1],
      condition = parts[2],
      avg_value = as.numeric(avg),
      pct_detected = as.numeric(pct),
      n_nuclei = length(cells),
      assay = assay_name,
      stringsAsFactors = FALSE
    )
  }) %>% bind_rows()

  out
}

rna_summary <- summarize_assay(obj, "RNA", slot_name = "data") %>%
  rename(avg_rna = avg_value, pct_rna_detected = pct_detected) %>%
  select(geneID, celltype, condition, avg_rna, pct_rna_detected, n_nuclei)

activity_summary <- summarize_assay(obj, "ACTIVITY", slot_name = "data") %>%
  rename(avg_activity = avg_value, pct_activity_detected = pct_detected) %>%
  select(geneID, celltype, condition, avg_activity, pct_activity_detected, n_nuclei)

write_tsv(rna_summary, file.path(OUT_DIR, "rna_by_celltype_condition.tsv"))
write_tsv(activity_summary, file.path(OUT_DIR, "activity_by_celltype_condition.tsv"))

# Calculate NaCl/Ctrl log2 fold changes where both conditions exist.
calc_log2fc <- function(df, value_col, pct_col, out_value_name) {
  df %>%
    select(geneID, celltype, condition, all_of(value_col), all_of(pct_col), n_nuclei) %>%
    pivot_wider(
      names_from = condition,
      values_from = c(all_of(value_col), all_of(pct_col), n_nuclei),
      values_fill = 0
    ) %>%
    mutate(
      log2FC_NaCl_vs_Ctrl = log2((.data[[paste0(value_col, "_Salinity")]] + 1e-6) /
                                   (.data[[paste0(value_col, "_Control")]] + 1e-6))
    ) %>%
    rename(!!out_value_name := log2FC_NaCl_vs_Ctrl)
}

# Support either "Salinity" or "NaCl" condition labels.
standardize_conditions <- function(df) {
  df %>%
    mutate(condition = case_when(
      condition %in% c("NaCl", "Salt", "salinity", "Salinity") ~ "Salinity",
      condition %in% c("Ctrl", "Ctl", "Control", "control") ~ "Control",
      TRUE ~ condition
    ))
}

rna_summary_std <- standardize_conditions(rna_summary)
activity_summary_std <- standardize_conditions(activity_summary)

rna_fc <- calc_log2fc(rna_summary_std, "avg_rna", "pct_rna_detected", "rna_log2FC_NaCl_vs_Ctrl")
activity_fc <- calc_log2fc(activity_summary_std, "avg_activity", "pct_activity_detected", "activity_log2FC_NaCl_vs_Ctrl")

write_tsv(rna_fc, file.path(OUT_DIR, "rna_log2fc_by_celltype.tsv"))
write_tsv(activity_fc, file.path(OUT_DIR, "activity_log2fc_by_celltype.tsv"))

message("Saved RNA and gene activity summaries by celltype × condition.")
