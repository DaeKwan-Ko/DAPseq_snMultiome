#!/usr/bin/env Rscript

# 02_parse_orthofinder_orthologues.R
# Parse OrthoFinder orthologue tables for sorghum ERF comparisons.
#
# This script is intentionally permissive because OrthoFinder filenames differ
# across runs. It searches ORTHOFINDER_DIR recursively for TSV files containing
# orthologue relationships and extracts rows involving sorghum ERF genes.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

if (!file.exists(ERF_CLADE_METADATA)) {
  stop("ERF clade metadata not found: ", ERF_CLADE_METADATA)
}

clade <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE)
if (!all(c("geneID", "clade") %in% colnames(clade))) {
  stop("ERF clade metadata must contain columns: geneID, clade")
}

orth_files <- list.files(
  ORTHOFINDER_DIR,
  pattern = "\\.tsv$|\\.txt$",
  recursive = TRUE,
  full.names = TRUE
)

if (length(orth_files) == 0) {
  message("No OrthoFinder TSV/TXT files found in: ", ORTHOFINDER_DIR)
  message("Skipping orthologue parsing.")
  quit(save = "no", status = 0)
}

parse_one_file <- function(path) {
  x <- tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) NULL
  )
  if (is.null(x) || ncol(x) < 2) return(tibble())

  # OrthoFinder Orthologues files commonly have columns:
  # Orthogroup, SpeciesA, SpeciesB
  # with comma-separated gene IDs.
  file_name <- basename(path)

  # Keep likely gene-containing columns.
  long <- x %>%
    mutate(row_id = row_number()) %>%
    pivot_longer(
      cols = -row_id,
      names_to = "source_column",
      values_to = "gene_list"
    ) %>%
    filter(!is.na(gene_list)) %>%
    separate_rows(gene_list, sep = ",\\s*|;\\s*") %>%
    mutate(
      geneID = str_trim(gene_list),
      file = file_name
    ) %>%
    filter(geneID != "")

  # Identify rows containing sorghum ERFs.
  sorghum_hits <- long %>%
    filter(str_detect(geneID, "Sobic\\.\\d{3}G\\d+")) %>%
    inner_join(clade, by = c("geneID" = "geneID")) %>%
    select(row_id, sorghum_geneID = geneID, clade, file)

  if (nrow(sorghum_hits) == 0) return(tibble())

  orthologs <- long %>%
    semi_join(sorghum_hits, by = "row_id") %>%
    filter(!str_detect(geneID, "Sobic\\.\\d{3}G\\d+")) %>%
    select(row_id, ortholog_geneID = geneID, ortholog_source_column = source_column, file)

  sorghum_hits %>%
    inner_join(orthologs, by = c("row_id", "file")) %>%
    mutate(
      comparison = case_when(
        str_detect(str_to_lower(file), "oryza|rice|osa|osativa") |
          str_detect(str_to_lower(ortholog_source_column), "oryza|rice|osa|osativa") ~ "sorghum_rice",
        str_detect(str_to_lower(file), "arabidopsis|athaliana|arabidopsis_thaliana|ath") |
          str_detect(str_to_lower(ortholog_source_column), "arabidopsis|athaliana|ath") ~ "sorghum_arabidopsis",
        str_detect(str_to_lower(file), "glycine|soybean|gmax|glycine_max") |
          str_detect(str_to_lower(ortholog_source_column), "glycine|soybean|gmax") ~ "sorghum_soybean",
        TRUE ~ "unknown"
      )
    ) %>%
    select(comparison, sorghum_geneID, ortholog_geneID, clade, source_file = file, ortholog_source_column)
}

ortholog_pairs <- map_dfr(orth_files, parse_one_file) %>%
  distinct() %>%
  arrange(comparison, clade, sorghum_geneID, ortholog_geneID)

if (nrow(ortholog_pairs) == 0) {
  message("No sorghum ERF orthologue pairs were parsed. Check OrthoFinder files and gene IDs.")
  quit(save = "no", status = 0)
}

out_file <- file.path(OUT_DIR, "ortholog_pairs_long.tsv")
write_tsv(ortholog_pairs, out_file)

summary_file <- file.path(OUT_DIR, "ortholog_pairs_summary.tsv")
ortholog_pairs %>%
  count(comparison, clade, name = "n_pairs") %>%
  write_tsv(summary_file)

message("Saved ortholog pairs: ", out_file)
message("Saved ortholog summary: ", summary_file)
message("Done.")
