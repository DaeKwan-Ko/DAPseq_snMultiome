#!/usr/bin/env Rscript

# 06_summarize_af3_metrics.R
# Summarize AlphaFold3 structural metrics across independent predictions.

source("scripts/05_alphafold3_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

designs <- NULL
if (file.exists(AF3_DESIGNS)) {
  designs <- read_tsv(AF3_DESIGNS, show_col_types = FALSE)
}

read_optional <- function(path) {
  if (file.exists(path)) {
    read_tsv(path, show_col_types = FALSE)
  } else {
    tibble()
  }
}

rmsd <- read_optional(DNA_RMSD)
contacts <- read_optional(CONTACTS)
beta2 <- read_optional(BETA2_DISTANCES)
conf <- read_optional(CONFIDENCE_SUMMARY)

# -------------------------------------------------------------------------
# DNA RMSD summary
# -------------------------------------------------------------------------

rmsd_summary <- tibble()
if (nrow(rmsd) > 0) {
  rmsd_summary <- rmsd %>%
    group_by(design_id) %>%
    summarise(
      n_models = n(),
      mean_dna_rmsd_A = mean(dna_rmsd_A, na.rm = TRUE),
      median_dna_rmsd_A = median(dna_rmsd_A, na.rm = TRUE),
      sd_dna_rmsd_A = sd(dna_rmsd_A, na.rm = TRUE),
      min_dna_rmsd_A = min(dna_rmsd_A, na.rm = TRUE),
      max_dna_rmsd_A = max(dna_rmsd_A, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------------------------------------------------------
# Contact summary
# -------------------------------------------------------------------------

contact_summary <- tibble()
if (nrow(contacts) > 0) {
  contact_summary <- contacts %>%
    group_by(design_id, motif_position) %>%
    summarise(
      n_models = n(),
      mean_base_atom_contacts = mean(n_protein_dna_base_atom_contacts, na.rm = TRUE),
      median_base_atom_contacts = median(n_protein_dna_base_atom_contacts, na.rm = TRUE),
      mean_contacting_residues = mean(n_contacting_protein_residues, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------------------------------------------------------
# Beta2 distance summary
# -------------------------------------------------------------------------

beta2_summary <- tibble()
if (nrow(beta2) > 0) {
  beta2_summary <- beta2 %>%
    group_by(design_id, beta2_site, motif_position) %>%
    summarise(
      n_models = n(),
      mean_min_distance_A = mean(min_beta2_to_dna_base_distance_A, na.rm = TRUE),
      median_min_distance_A = median(min_beta2_to_dna_base_distance_A, na.rm = TRUE),
      min_distance_A = min(min_beta2_to_dna_base_distance_A, na.rm = TRUE),
      max_distance_A = max(min_beta2_to_dna_base_distance_A, na.rm = TRUE),
      direct_contact_fraction_4A = mean(min_beta2_to_dna_base_distance_A <= 4, na.rm = TRUE),
      .groups = "drop"
    )
}

# -------------------------------------------------------------------------
# Combined source-data table
# -------------------------------------------------------------------------

combined <- rmsd_summary

if (nrow(designs) > 0 && nrow(combined) > 0) {
  combined <- combined %>% left_join(designs, by = "design_id")
}

write_tsv(rmsd_summary, file.path(OUT_DIR, "af3_dna_rmsd_summary.tsv"))
write_tsv(contact_summary, file.path(OUT_DIR, "af3_base_contact_summary.tsv"))
write_tsv(beta2_summary, file.path(OUT_DIR, "af3_beta2_distance_summary.tsv"))

# Save long source data for Extended Data Fig. 4-style panels.
if (nrow(rmsd) > 0) {
  rmsd_source <- rmsd %>%
    mutate(metric = "DNA_RMSD_A", value = dna_rmsd_A) %>%
    select(design_id, model_id, metric, value, everything())
  write_tsv(rmsd_source, file.path(OUT_DIR, "af3_extended_data_fig4_dna_rmsd_source_data.tsv"))
}

if (nrow(contacts) > 0) {
  contact_source <- contacts %>%
    mutate(metric = "protein_DNA_base_atom_contacts", value = n_protein_dna_base_atom_contacts) %>%
    select(design_id, model_id, motif_position, metric, value, everything())
  write_tsv(contact_source, file.path(OUT_DIR, "af3_extended_data_fig4_contact_source_data.tsv"))
}

if (nrow(beta2) > 0) {
  beta2_source <- beta2 %>%
    mutate(metric = "beta2_to_DNA_base_min_distance_A", value = min_beta2_to_dna_base_distance_A) %>%
    select(design_id, model_id, beta2_site, motif_position, metric, value, everything())
  write_tsv(beta2_source, file.path(OUT_DIR, "af3_extended_data_fig4_beta2_distance_source_data.tsv"))
}

# Compact overview
overview <- list(
  rmsd_summary = rmsd_summary,
  contact_summary = contact_summary,
  beta2_summary = beta2_summary
)

message("Saved AlphaFold3 metric summaries to: ", OUT_DIR)
message("RMSD rows: ", nrow(rmsd))
message("Contact rows: ", nrow(contacts))
message("Beta2 distance rows: ", nrow(beta2))
message("Done.")
