#!/usr/bin/env Rscript

# 04_plot_network_summaries.R
# Plot Cortex_1 ERF network summaries and driver rankings.

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
# Weighted out-degree by ERF
# -------------------------------------------------------------------------

outdegree <- read_optional_tsv(WEIGHTED_OUT_DEGREE)

if (nrow(outdegree) > 0 && all(c("erf_geneID", "weighted_out_degree") %in% colnames(outdegree))) {
  outdegree_plot <- outdegree %>%
    mutate(
      clade_full = if ("erf_clade_full" %in% colnames(.)) erf_clade_full else standardize_clade(clade),
      erf_geneID = fct_reorder(erf_geneID, weighted_out_degree)
    )

  write_tsv(outdegree_plot, file.path(SOURCE_DATA_DIR, "fig4_weighted_out_degree_cortex1_source_data.tsv"))

  p <- ggplot(outdegree_plot, aes(x = erf_geneID, y = weighted_out_degree)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_grid(clade_full ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "ERF",
      y = "Weighted out-degree",
      title = paste0("ERF network connectivity in ", TARGET_CELLTYPE)
    ) +
    theme_manuscript(base_size = 7)

  save_plot(p, file.path(PLOT_DIR, "fig4_weighted_out_degree_cortex1.pdf"), width = 4.2, height = 5.2)
}

# -------------------------------------------------------------------------
# Driver ranking
# -------------------------------------------------------------------------

ranking <- read_optional_tsv(DRIVER_RANKING)

if (nrow(ranking) > 0 && all(c("erf_geneID", "driver_score", "clade") %in% colnames(ranking))) {
  ranking_plot <- ranking %>%
    group_by(clade) %>%
    arrange(clade_rank, .by_group = TRUE) %>%
    slice_head(n = 10) %>%
    ungroup() %>%
    mutate(
      clade_full = if ("erf_clade_full" %in% colnames(.)) erf_clade_full else standardize_clade(clade),
      erf_geneID = fct_reorder(erf_geneID, driver_score)
    )

  write_tsv(ranking_plot, file.path(SOURCE_DATA_DIR, "fig4_erf_driver_ranking_cortex1_source_data.tsv"))

  p <- ggplot(ranking_plot, aes(x = erf_geneID, y = driver_score)) +
    geom_col(width = 0.75) +
    coord_flip() +
    facet_grid(clade_full ~ ., scales = "free_y", space = "free_y") +
    labs(
      x = "ERF",
      y = "Matched-clade driver score",
      title = paste0("Candidate ERF driver ranking in ", TARGET_CELLTYPE)
    ) +
    theme_manuscript(base_size = 7)

  save_plot(p, file.path(PLOT_DIR, "fig4_erf_driver_ranking_cortex1.pdf"), width = 4.0, height = 4.2)
}

# -------------------------------------------------------------------------
# Edge weight distribution
# -------------------------------------------------------------------------

edges <- read_optional_tsv(NETWORK_EDGES)

if (nrow(edges) > 0 && "edge_weight" %in% colnames(edges)) {
  edge_plot <- edges %>%
    mutate(
      clade_full = if ("erf_clade_full" %in% colnames(.)) erf_clade_full else standardize_clade(clade)
    )

  write_tsv(edge_plot, file.path(SOURCE_DATA_DIR, "fig4_network_edges_cortex1_source_data.tsv"))

  p <- ggplot(edge_plot, aes(x = clade_full, y = edge_weight)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 0.6, alpha = 0.45) +
    labs(
      x = "ERF clade",
      y = "Edge weight",
      title = paste0("ERF-to-target edge weights in ", TARGET_CELLTYPE)
    ) +
    theme_manuscript()

  save_plot(p, file.path(PLOT_DIR, "fig4_network_edge_weight_distribution.pdf"), width = 2.8, height = 3.0)
}

message("Finished network summary plotting.")
