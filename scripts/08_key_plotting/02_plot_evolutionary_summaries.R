#!/usr/bin/env Rscript

# 02_plot_evolutionary_summaries.R
# Plot Ka/Ks and cross-lineage beta2 signature summaries.

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
# Pooled Ka/Ks boxplot
# -------------------------------------------------------------------------

kaks <- read_optional_tsv(KAKS_PAIRS)

if (nrow(kaks) > 0 && all(c("clade", "kaks") %in% colnames(kaks))) {
  kaks_plot <- kaks %>%
    mutate(
      clade_full = if ("clade_full" %in% colnames(.)) clade_full else standardize_clade(clade),
      comparison = ifelse(is.na(comparison), "unknown", comparison)
    ) %>%
    filter(is.finite(kaks))

  write_tsv(kaks_plot, file.path(SOURCE_DATA_DIR, "fig3_kaks_pooled_boxplot_source_data.tsv"))

  p <- ggplot(kaks_plot, aes(x = clade_full, y = kaks)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 0.7, alpha = 0.5) +
    labs(
      x = "ERF clade",
      y = "Ka/Ks across AP2/ERF DBD",
      title = "Pooled AP2/ERF DBD Ka/Ks"
    ) +
    theme_manuscript()

  save_plot(p, file.path(PLOT_DIR, "fig3_kaks_pooled_boxplot.pdf"), width = 2.4, height = 3.0)

  p2 <- ggplot(kaks_plot, aes(x = clade_full, y = kaks)) +
    geom_boxplot(outlier.shape = NA, width = 0.55) +
    geom_jitter(width = 0.12, size = 0.6, alpha = 0.45) +
    facet_wrap(~comparison, scales = "free_y") +
    labs(
      x = "ERF clade",
      y = "Ka/Ks",
      title = "AP2/ERF DBD Ka/Ks by ortholog comparison"
    ) +
    theme_manuscript()

  save_plot(p2, file.path(PLOT_DIR, "fig3_kaks_by_comparison_boxplot.pdf"), width = 5.0, height = 3.0)
}

# -------------------------------------------------------------------------
# Cross-lineage beta2 stacked bar
# -------------------------------------------------------------------------

beta2 <- read_optional_tsv(BETA2_LINEAGE_SUMMARY)

if (nrow(beta2) > 0 && all(c("species", "lineage_order", "n_AA", "n_WV", "n_neither") %in% colnames(beta2))) {
  beta2_long <- beta2 %>%
    select(lineage_order, species, taxonomy, major_group, n_AA, n_WV, n_neither) %>%
    pivot_longer(
      cols = c(n_AA, n_WV, n_neither),
      names_to = "signature_class",
      values_to = "n_erfs"
    ) %>%
    mutate(
      signature_class = recode(signature_class, n_AA = "AA", n_WV = "WV", n_neither = "neither"),
      species = fct_reorder(species, lineage_order)
    )

  write_tsv(beta2_long, file.path(SOURCE_DATA_DIR, "fig3_cross_lineage_beta2_stacked_bar_source_data.tsv"))

  p <- ggplot(beta2_long, aes(x = species, y = n_erfs, fill = signature_class)) +
    geom_col(width = 0.75) +
    coord_flip() +
    labs(
      x = NULL,
      y = "Number of ERFs",
      fill = "β2 signature",
      title = "Cross-lineage distribution of ERF β2 signatures"
    ) +
    theme_manuscript()

  save_plot(p, file.path(PLOT_DIR, "fig3_cross_lineage_beta2_stacked_bar.pdf"), width = 4.4, height = 3.6)
}

message("Finished evolutionary summary plotting.")
