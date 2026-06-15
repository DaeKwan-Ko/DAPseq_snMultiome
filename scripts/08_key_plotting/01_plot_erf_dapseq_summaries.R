#!/usr/bin/env Rscript

# 01_plot_erf_dapseq_summaries.R
# Plot ERF DAP-seq summary panels, including clade-biased target counts and
# peak annotation distributions when available.

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
# ERF clade-biased recurrent target counts
# -------------------------------------------------------------------------

clade_counts <- read_optional_tsv(ERF_TARGET_CLADE_COUNTS)

if (nrow(clade_counts) > 0) {
  if ("assigned_geneID" %in% colnames(clade_counts) && !"geneID" %in% colnames(clade_counts)) {
    clade_counts <- clade_counts %>% rename(geneID = assigned_geneID)
  }

  if ("target_class" %in% colnames(clade_counts)) {
    target_class_counts <- clade_counts %>%
      filter(target_class %in% c("ECA_biased_recurrent", "ECB_biased_recurrent", "other")) %>%
      count(target_class, name = "n_targets") %>%
      mutate(
        target_class = factor(
          target_class,
          levels = c("ECA_biased_recurrent", "ECB_biased_recurrent", "other")
        )
      )

    write_tsv(target_class_counts, file.path(SOURCE_DATA_DIR, "fig2_erf_clade_target_counts_source_data.tsv"))

    p <- ggplot(target_class_counts, aes(x = target_class, y = n_targets)) +
      geom_col(width = 0.7) +
      geom_text(aes(label = n_targets), vjust = -0.3, size = 2.5) +
      labs(
        x = NULL,
        y = "Number of candidate target genes",
        title = "ERF clade-biased recurrent targets"
      ) +
      theme_manuscript() +
      theme(axis.text.x = element_text(angle = 35, hjust = 1))

    save_plot(p, file.path(PLOT_DIR, "fig2_erf_clade_target_counts.pdf"), width = 3.6, height = 3.0)
  }

  if (all(c("n_clade_A", "n_clade_B") %in% colnames(clade_counts))) {
    scatter_data <- clade_counts %>%
      mutate(
        target_class = ifelse(is.na(target_class), "other", target_class)
      )

    write_tsv(scatter_data, file.path(SOURCE_DATA_DIR, "fig2_erf_target_clade_count_scatter_source_data.tsv"))

    p <- ggplot(scatter_data, aes(x = n_clade_A, y = n_clade_B)) +
      geom_point(aes(shape = target_class), alpha = 0.7, size = 1.5) +
      labs(
        x = "Number of ECA ERFs associated with target",
        y = "Number of ECB ERFs associated with target",
        title = "ERF clade association of candidate targets"
      ) +
      theme_manuscript()

    save_plot(p, file.path(PLOT_DIR, "fig2_erf_target_clade_count_scatter.pdf"), width = 3.4, height = 3.2)
  }
}

# -------------------------------------------------------------------------
# Peak annotation distributions by TF family
# -------------------------------------------------------------------------

peak_summary <- read_optional_tsv(PEAK_ANNOT_SUMMARY)

if (nrow(peak_summary) > 0 && all(c("tf_geneID", "annotation_class", "percent_peaks") %in% colnames(peak_summary))) {
  peak_plot_data <- peak_summary %>%
    mutate(
      tf_family = ifelse("tf_family" %in% colnames(.), tf_family, NA_character_)
    )

  write_tsv(peak_plot_data, file.path(SOURCE_DATA_DIR, "peak_annotation_summary_source_data.tsv"))

  p <- ggplot(peak_plot_data, aes(x = tf_geneID, y = percent_peaks, fill = annotation_class)) +
    geom_col(width = 0.85) +
    facet_grid(. ~ tf_family, scales = "free_x", space = "free_x") +
    labs(
      x = "TF",
      y = "Peaks (%)",
      fill = "Annotation",
      title = "Genomic annotation of DAP-seq peaks"
    ) +
    theme_manuscript(base_size = 7) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.position = "right"
    )

  save_plot(p, file.path(PLOT_DIR, "extended_peak_annotation_by_tf_family.pdf"), width = 7.0, height = 3.0)
}

message("Finished ERF/DAP-seq summary plotting.")
