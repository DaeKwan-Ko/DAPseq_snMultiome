#!/usr/bin/env Rscript

# 06_build_erf_target_matrices.R
# Build ERF target matrices and identify clade-biased recurrent targets.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

target_file <- file.path(OUT_DIR, "dapseq_candidate_targets_long.tsv")
if (!file.exists(target_file)) {
  stop("Candidate target file not found. Run 04_define_candidate_targets.R first: ", target_file)
}

if (!file.exists(ERF_CLADE_METADATA)) {
  stop("ERF_CLADE_METADATA not found: ", ERF_CLADE_METADATA)
}

targets <- read_tsv(target_file, show_col_types = FALSE)
clade <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE)

if (!all(c("geneID", "clade") %in% colnames(clade))) {
  stop("ERF clade metadata must contain columns: geneID, clade")
}

erf_targets <- targets %>%
  inner_join(clade, by = c("tf_geneID" = "geneID")) %>%
  distinct(tf_geneID, clade, assigned_geneID)

target_matrix <- erf_targets %>%
  mutate(bound = 1L) %>%
  select(assigned_geneID, tf_geneID, bound) %>%
  pivot_wider(
    names_from = tf_geneID,
    values_from = bound,
    values_fill = 0
  ) %>%
  arrange(assigned_geneID)

clade_counts <- erf_targets %>%
  count(assigned_geneID, clade, name = "n_erfs") %>%
  pivot_wider(
    names_from = clade,
    values_from = n_erfs,
    values_fill = 0,
    names_prefix = "n_clade_"
  ) %>%
  mutate(
    n_clade_A = ifelse("n_clade_A" %in% colnames(.), n_clade_A, 0),
    n_clade_B = ifelse("n_clade_B" %in% colnames(.), n_clade_B, 0),
    n_total_erfs = n_clade_A + n_clade_B,
    target_class = case_when(
      n_clade_A >= 6 & n_clade_B <= 3 ~ "ECA_biased_recurrent",
      n_clade_B >= 6 & n_clade_A <= 3 ~ "ECB_biased_recurrent",
      TRUE ~ "other"
    )
  ) %>%
  arrange(desc(n_total_erfs), assigned_geneID)

biased_recurrent <- clade_counts %>%
  filter(target_class %in% c("ECA_biased_recurrent", "ECB_biased_recurrent"))

write_tsv(erf_targets, file.path(OUT_DIR, "erf_candidate_targets_long.tsv"))
write_tsv(target_matrix, file.path(OUT_DIR, "erf_target_matrix.tsv"))
write_tsv(clade_counts, file.path(OUT_DIR, "erf_target_clade_counts.tsv"))
write_tsv(biased_recurrent, file.path(OUT_DIR, "erf_clade_biased_recurrent_targets.tsv"))

message("Saved ERF target matrices and clade-biased recurrent target table.")
message("Number of ERF target genes: ", n_distinct(erf_targets$assigned_geneID))
message("Number of ECA-biased recurrent targets: ", sum(biased_recurrent$target_class == "ECA_biased_recurrent"))
message("Number of ECB-biased recurrent targets: ", sum(biased_recurrent$target_class == "ECB_biased_recurrent"))
message("Done.")
