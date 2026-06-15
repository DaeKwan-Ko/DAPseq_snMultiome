#!/usr/bin/env Rscript

# 03_summarize_peak_annotations.R
# Summarize ChIPseeker annotation categories for each TF.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

annotation_file <- file.path(OUT_DIR, "all_peak_annotations.tsv")
if (!file.exists(annotation_file)) {
  stop("Annotation file not found. Run 02_annotate_peaks_chipseeker.R first: ", annotation_file)
}

anno <- read_tsv(annotation_file, show_col_types = FALSE)

tf_meta <- NULL
if (file.exists(TF_METADATA)) {
  tf_meta <- read_tsv(TF_METADATA, show_col_types = FALSE)
  # Accept either geneID/tf_family or older label/tf headers
  if (all(c("label", "tf") %in% colnames(tf_meta)) && !("geneID" %in% colnames(tf_meta))) {
    tf_meta <- tf_meta %>% rename(geneID = label, tf_family = tf)
  }
}

summary_by_tf <- anno %>%
  count(tf_geneID, annotation_class, name = "n_peaks") %>%
  group_by(tf_geneID) %>%
  mutate(
    total_peaks = sum(n_peaks),
    percent_peaks = 100 * n_peaks / total_peaks
  ) %>%
  ungroup()

if (!is.null(tf_meta) && all(c("geneID", "tf_family") %in% colnames(tf_meta))) {
  summary_by_tf <- summary_by_tf %>%
    left_join(tf_meta %>% select(geneID, tf_family), by = c("tf_geneID" = "geneID")) %>%
    relocate(tf_family, .after = tf_geneID)
}

summary_wide <- summary_by_tf %>%
  select(tf_geneID, tf_family, annotation_class, percent_peaks) %>%
  pivot_wider(
    names_from = annotation_class,
    values_from = percent_peaks,
    values_fill = 0,
    names_prefix = "percent_"
  )

write_tsv(summary_by_tf, file.path(OUT_DIR, "peak_annotation_summary_by_tf.tsv"))
write_tsv(summary_wide, file.path(OUT_DIR, "peak_annotation_summary_by_tf_wide.tsv"))

# Family-level summary
if ("tf_family" %in% colnames(summary_by_tf)) {
  family_summary <- summary_by_tf %>%
    group_by(tf_family, annotation_class) %>%
    summarise(
      mean_percent_peaks = mean(percent_peaks, na.rm = TRUE),
      median_percent_peaks = median(percent_peaks, na.rm = TRUE),
      n_tfs = n_distinct(tf_geneID),
      .groups = "drop"
    )

  write_tsv(family_summary, file.path(OUT_DIR, "peak_annotation_summary_by_family.tsv"))
}

message("Saved peak annotation summaries to: ", OUT_DIR)
message("Done.")
