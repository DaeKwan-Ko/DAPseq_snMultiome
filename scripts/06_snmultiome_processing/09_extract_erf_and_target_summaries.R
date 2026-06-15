#!/usr/bin/env Rscript

# 09_extract_erf_and_target_summaries.R
# Extract ERF and candidate target gene RNA/activity summaries.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

rna_file <- file.path(OUT_DIR, "rna_by_celltype_condition.tsv")
activity_file <- file.path(OUT_DIR, "activity_by_celltype_condition.tsv")
rna_fc_file <- file.path(OUT_DIR, "rna_log2fc_by_celltype.tsv")
activity_fc_file <- file.path(OUT_DIR, "activity_log2fc_by_celltype.tsv")

if (!file.exists(rna_file) || !file.exists(activity_file)) {
  stop("RNA/activity summary files not found. Run 08_summarize_rna_activity_by_group.R first.")
}

rna <- read_tsv(rna_file, show_col_types = FALSE)
activity <- read_tsv(activity_file, show_col_types = FALSE)
rna_fc <- if (file.exists(rna_fc_file)) read_tsv(rna_fc_file, show_col_types = FALSE) else tibble()
activity_fc <- if (file.exists(activity_fc_file)) read_tsv(activity_fc_file, show_col_types = FALSE) else tibble()

# ERF metadata
erf_meta <- tibble()
if (file.exists(ERF_CLADE_METADATA)) {
  erf_meta <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE) %>%
    mutate(clade_full = recode(clade, A = "ECA", B = "ECB", .default = clade))
}

salinity <- tibble()
if (file.exists(SALINITY_ERF_METADATA)) {
  salinity <- read_tsv(SALINITY_ERF_METADATA, show_col_types = FALSE) %>%
    mutate(salinity_linked = TRUE) %>%
    select(geneID, salinity_linked)
}

if (nrow(erf_meta) > 0) {
  erf_summary <- rna %>%
    inner_join(erf_meta, by = "geneID") %>%
    left_join(activity, by = c("geneID", "celltype", "condition", "n_nuclei")) %>%
    left_join(salinity, by = "geneID") %>%
    mutate(salinity_linked = ifelse(is.na(salinity_linked), FALSE, salinity_linked)) %>%
    arrange(clade, geneID, celltype, condition)

  write_tsv(erf_summary, file.path(OUT_DIR, "erf_rna_activity_summary_by_celltype_condition.tsv"))

  if (nrow(rna_fc) > 0 && nrow(activity_fc) > 0) {
    erf_fc <- rna_fc %>%
      inner_join(erf_meta, by = "geneID") %>%
      left_join(activity_fc %>% select(geneID, celltype, activity_log2FC_NaCl_vs_Ctrl), by = c("geneID", "celltype")) %>%
      left_join(salinity, by = "geneID") %>%
      mutate(salinity_linked = ifelse(is.na(salinity_linked), FALSE, salinity_linked)) %>%
      arrange(clade, geneID, celltype)

    write_tsv(erf_fc, file.path(OUT_DIR, "erf_rna_activity_log2fc_by_celltype.tsv"))
  }
}

# Candidate target genes
if (file.exists(CANDIDATE_TARGET_GENES)) {
  targets <- read_tsv(CANDIDATE_TARGET_GENES, show_col_types = FALSE)
  if (!"geneID" %in% colnames(targets)) {
    stop("candidate_target_genes.tsv must contain a geneID column")
  }

  target_summary <- rna %>%
    inner_join(targets, by = "geneID") %>%
    left_join(activity, by = c("geneID", "celltype", "condition", "n_nuclei")) %>%
    arrange(geneID, celltype, condition)

  write_tsv(target_summary, file.path(OUT_DIR, "target_gene_rna_activity_summary_by_celltype_condition.tsv"))

  if (nrow(rna_fc) > 0 && nrow(activity_fc) > 0) {
    target_fc <- rna_fc %>%
      inner_join(targets, by = "geneID") %>%
      left_join(activity_fc %>% select(geneID, celltype, activity_log2FC_NaCl_vs_Ctrl), by = c("geneID", "celltype")) %>%
      arrange(geneID, celltype)

    write_tsv(target_fc, file.path(OUT_DIR, "target_gene_rna_activity_log2fc_by_celltype.tsv"))
  }
} else {
  message("Candidate target metadata not found: ", CANDIDATE_TARGET_GENES)
  message("Skipping target gene summary.")
}

message("Saved ERF and candidate target RNA/activity summaries.")
