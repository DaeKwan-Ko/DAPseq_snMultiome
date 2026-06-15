#!/usr/bin/env Rscript

# 05_build_cross_species_beta2_summary.R
# Summarize cross-species ERF β2 signatures across representative green plant lineages.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

species_file <- file.path(OUT_DIR, "species_metadata_prepared.tsv")
if (!file.exists(species_file)) {
  stop("Prepared species metadata not found. Run 01_prepare_species_metadata.R first.")
}

if (!file.exists(CROSS_SPECIES_BETA2)) {
  message("CROSS_SPECIES_BETA2 not found: ", CROSS_SPECIES_BETA2)
  message("Expected columns: species, geneID, beta2_signature, beta2_signature_class")
  message("Skipping cross-species beta2 summary.")
  quit(save = "no", status = 0)
}

species_meta <- read_tsv(species_file, show_col_types = FALSE)
beta2 <- read_tsv(CROSS_SPECIES_BETA2, show_col_types = FALSE)

required <- c("species", "geneID", "beta2_signature_class")
missing <- setdiff(required, colnames(beta2))
if (length(missing) > 0) {
  stop("Missing required columns in CROSS_SPECIES_BETA2: ", paste(missing, collapse = ", "))
}

beta2 <- beta2 %>%
  mutate(
    beta2_signature_class = case_when(
      beta2_signature_class %in% c("AA", "WV", "neither") ~ beta2_signature_class,
      beta2_signature == "AA" ~ "AA",
      beta2_signature == "WV" ~ "WV",
      TRUE ~ "neither"
    ),
    species = str_replace_all(species, "_", " ")
  )

# Normalize some common species-name variants.
normalize_species <- function(x) {
  x %>%
    str_replace_all("Oryza Sativa", "Oryza sativa") %>%
    str_replace_all("Glycine Max", "Glycine max") %>%
    str_replace_all("Arabidopsis Thaliana", "Arabidopsis thaliana")
}

beta2 <- beta2 %>% mutate(species = normalize_species(species))
species_meta <- species_meta %>% mutate(species = normalize_species(species))

counts <- beta2 %>%
  count(species, beta2_signature_class, name = "n_erfs") %>%
  complete(
    species,
    beta2_signature_class = c("AA", "WV", "neither"),
    fill = list(n_erfs = 0)
  ) %>%
  left_join(species_meta, by = "species") %>%
  arrange(lineage_order, beta2_signature_class)

summary_wide <- counts %>%
  select(lineage_order, species, taxonomy, major_group, total_erfs_analyzed, beta2_signature_class, n_erfs) %>%
  pivot_wider(
    names_from = beta2_signature_class,
    values_from = n_erfs,
    values_fill = 0,
    names_prefix = "n_"
  ) %>%
  mutate(
    n_total_classified = n_AA + n_WV + n_neither,
    pct_AA = 100 * n_AA / n_total_classified,
    pct_WV = 100 * n_WV / n_total_classified,
    pct_neither = 100 * n_neither / n_total_classified,
    wv_detected = n_WV > 0,
    aa_detected = n_AA > 0
  ) %>%
  arrange(lineage_order)

write_tsv(beta2, file.path(OUT_DIR, "cross_species_beta2_signature_long.tsv"))
write_tsv(counts, file.path(OUT_DIR, "cross_species_beta2_signature_counts.tsv"))
write_tsv(summary_wide, file.path(OUT_DIR, "cross_lineage_beta2_summary_for_plotting.tsv"))

message("Saved cross-species beta2 summaries.")
message("Total ERFs classified: ", nrow(beta2))
print(summary_wide %>% select(species, n_AA, n_WV, n_neither, n_total_classified))
message("Done.")
