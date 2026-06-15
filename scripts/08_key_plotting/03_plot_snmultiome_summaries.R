#!/usr/bin/env Rscript

# 03_plot_snmultiome_summaries.R
# Plot snMultiome celltype composition and RNA/activity summaries.

source("scripts/08_key_plotting/00_config.R")
source("scripts/08_key_plotting/plotting_helpers.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
})

# -------------------------------------------------------------------------
# Celltype × condition nucleus counts
# -------------------------------------------------------------------------

counts <- read_optional_tsv(CELLTYPE_COUNTS)

if (nrow(counts) > 0 && all(c("celltype", "condition", "n_nuclei") %in% colnames(counts))) {
  counts_plot <- counts %>%
    mutate(
      condition = standardize_condition(condition),
      celltype = fct_reorder(celltype, n_nuclei, .fun = sum)
    )

  write_tsv(counts_plot, file.path(SOURCE_DATA_DIR, "fig4_celltype_condition_nucleus_counts_source_data.tsv"))

  p <- ggplot(counts_plot, aes(x = celltype, y = n_nuclei, fill = condition)) +
    geom_col(position = position_dodge(width = 0.75), width = 0.7) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Number of nuclei",
      fill = "Condition",
      title = "snMultiome nuclei by cell type and condition"
    ) +
    theme_manuscript()

  save_plot(p, file.path(PLOT_DIR, "fig4_celltype_condition_nucleus_counts.pdf"), width = 3.6, height = 3.2)
}

# -------------------------------------------------------------------------
# ERF RNA log2FC by cell type and clade
# -------------------------------------------------------------------------

erf_fc <- read_optional_tsv(ERF_RNA_ACTIVITY_FC)

if (nrow(erf_fc) > 0 && all(c("geneID", "celltype", "clade", "rna_log2FC_NaCl_vs_Ctrl") %in% colnames(erf_fc))) {
  erf_plot <- erf_fc %>%
    mutate(
      clade_full = if ("clade_full" %in% colnames(.)) clade_full else standardize_clade(clade),
      salinity_linked = if ("salinity_linked" %in% colnames(.)) salinity_linked else NA
    )

  write_tsv(erf_plot, file.path(SOURCE_DATA_DIR, "fig4_erf_rna_log2fc_by_celltype_source_data.tsv"))

  p <- ggplot(erf_plot, aes(x = clade_full, y = rna_log2FC_NaCl_vs_Ctrl)) +
    geom_violin(width = 0.8, trim = FALSE) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.5) +
    facet_wrap(~celltype, scales = "free_y") +
    labs(
      x = "ERF clade",
      y = "RNA log2FC (NaCl/Ctrl)",
      title = "Cell type-resolved salinity response of ERFs"
    ) +
    theme_manuscript()

  save_plot(p, file.path(PLOT_DIR, "fig4_erf_rna_log2fc_by_celltype.pdf"), width = 6.0, height = 3.4)
}

# -------------------------------------------------------------------------
# Target RNA/activity log2FC in Cortex_1
# -------------------------------------------------------------------------

target_fc <- read_optional_tsv(TARGET_RNA_ACTIVITY_FC)

if (nrow(target_fc) > 0 && "celltype" %in% colnames(target_fc)) {
  value_cols <- intersect(
    c("rna_log2FC_NaCl_vs_Ctrl", "activity_log2FC_NaCl_vs_Ctrl"),
    colnames(target_fc)
  )

  if (length(value_cols) > 0) {
    target_plot <- target_fc %>%
      filter(celltype == TARGET_CELLTYPE) %>%
      mutate(target_group = if ("target_group" %in% colnames(.)) target_group else "candidate_target") %>%
      select(geneID, celltype, target_group, all_of(value_cols)) %>%
      pivot_longer(
        cols = all_of(value_cols),
        names_to = "assay",
        values_to = "log2FC_NaCl_vs_Ctrl"
      ) %>%
      mutate(
        assay = recode(
          assay,
          rna_log2FC_NaCl_vs_Ctrl = "RNA",
          activity_log2FC_NaCl_vs_Ctrl = "ATAC-derived gene activity"
        )
      )

    write_tsv(target_plot, file.path(SOURCE_DATA_DIR, "fig4_target_rna_activity_cortex1_source_data.tsv"))

    p <- ggplot(target_plot, aes(x = target_group, y = log2FC_NaCl_vs_Ctrl)) +
      geom_violin(width = 0.8, trim = FALSE) +
      geom_jitter(width = 0.12, size = 0.55, alpha = 0.45) +
      facet_wrap(~assay, scales = "free_y") +
      labs(
        x = NULL,
        y = "log2FC (NaCl/Ctrl)",
        title = paste0("Clade-biased recurrent targets in ", TARGET_CELLTYPE)
      ) +
      theme_manuscript() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))

    save_plot(p, file.path(PLOT_DIR, "fig4_target_rna_activity_cortex1.pdf"), width = 5.0, height = 3.2)
  }
}

message("Finished snMultiome summary plotting.")
