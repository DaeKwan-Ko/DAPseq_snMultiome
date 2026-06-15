#!/usr/bin/env Rscript

# 04_summarize_kaks_by_clade.R
# Summarize AP2/ERF DBD Ka/Ks values by ERF clade and species comparison.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
})

annotated_file <- file.path(OUT_DIR, "kaks_ortholog_pairs_annotated.tsv")
if (!file.exists(annotated_file)) {
  message("Annotated Ka/Ks file not found. Run 03_prepare_kaks_ortholog_pairs.R first.")
  quit(save = "no", status = 0)
}

kaks <- read_tsv(annotated_file, show_col_types = FALSE) %>%
  filter(is.finite(kaks), !is.na(clade))

summary_by_clade <- kaks %>%
  group_by(comparison, clade, clade_full) %>%
  summarise(
    n_pairs = n(),
    n_sorghum_genes = n_distinct(sorghum_geneID),
    median_kaks = median(kaks, na.rm = TRUE),
    mean_kaks = mean(kaks, na.rm = TRUE),
    q1_kaks = quantile(kaks, 0.25, na.rm = TRUE),
    q3_kaks = quantile(kaks, 0.75, na.rm = TRUE),
    min_kaks = min(kaks, na.rm = TRUE),
    max_kaks = max(kaks, na.rm = TRUE),
    .groups = "drop"
  )

summary_pooled <- kaks %>%
  mutate(comparison = "pooled") %>%
  group_by(comparison, clade, clade_full) %>%
  summarise(
    n_pairs = n(),
    n_sorghum_genes = n_distinct(sorghum_geneID),
    median_kaks = median(kaks, na.rm = TRUE),
    mean_kaks = mean(kaks, na.rm = TRUE),
    q1_kaks = quantile(kaks, 0.25, na.rm = TRUE),
    q3_kaks = quantile(kaks, 0.75, na.rm = TRUE),
    min_kaks = min(kaks, na.rm = TRUE),
    max_kaks = max(kaks, na.rm = TRUE),
    .groups = "drop"
  )

summary_all <- bind_rows(summary_by_clade, summary_pooled)

# Wilcoxon tests: ECA vs ECB for each comparison and pooled
run_wilcox <- function(df) {
  if (n_distinct(df$clade) < 2) {
    return(tibble(p_value = NA_real_, statistic = NA_real_, n_A = sum(df$clade == "A"), n_B = sum(df$clade == "B")))
  }
  test <- wilcox.test(kaks ~ clade, data = df, exact = FALSE)
  tibble(
    p_value = test$p.value,
    statistic = unname(test$statistic),
    n_A = sum(df$clade == "A"),
    n_B = sum(df$clade == "B")
  )
}

tests_by_comparison <- kaks %>%
  group_by(comparison) %>%
  group_modify(~run_wilcox(.x)) %>%
  ungroup()

tests_pooled <- kaks %>%
  mutate(comparison = "pooled") %>%
  group_by(comparison) %>%
  group_modify(~run_wilcox(.x)) %>%
  ungroup()

tests_all <- bind_rows(tests_by_comparison, tests_pooled) %>%
  mutate(p_adjust_bh = p.adjust(p_value, method = "BH"))

write_tsv(summary_all, file.path(OUT_DIR, "kaks_summary_by_clade.tsv"))
write_tsv(tests_all, file.path(OUT_DIR, "kaks_wilcoxon_tests.tsv"))

# Save a compact source-data table for box plots.
plot_data <- kaks %>%
  select(comparison, sorghum_geneID, ortholog_geneID, clade, clade_full, ka, ks, kaks)

write_tsv(plot_data, file.path(OUT_DIR, "kaks_boxplot_source_data.tsv"))

message("Saved Ka/Ks summaries and Wilcoxon test tables.")
message("Pooled usable pairs: ", nrow(kaks))
message("Done.")
