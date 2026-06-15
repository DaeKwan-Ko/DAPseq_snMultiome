#!/usr/bin/env Rscript

# 05_go_enrichment_hypergeometric.R
# Perform simple hypergeometric GO enrichment for DAP-seq-associated candidate targets.
#
# Required input:
#   metadata/sorghum_gene2go.tsv
# Required columns:
#   geneID, go_id
# Optional columns:
#   go_name, ontology
#
# Output:
#   go_enrichment_by_tf.tsv
#   go_enrichment_erf_clade_targets.tsv

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(purrr)
})

target_file <- file.path(OUT_DIR, "dapseq_candidate_targets_long.tsv")
if (!file.exists(target_file)) {
  stop("Candidate target file not found. Run 04_define_candidate_targets.R first: ", target_file)
}

if (!file.exists(GENE2GO_FILE)) {
  message("GENE2GO_FILE not found: ", GENE2GO_FILE)
  message("Skipping GO enrichment. Create metadata/sorghum_gene2go.tsv and rerun this script.")
  quit(save = "no", status = 0)
}

targets <- read_tsv(target_file, show_col_types = FALSE)
gene2go <- read_tsv(GENE2GO_FILE, show_col_types = FALSE)

if (!all(c("geneID", "go_id") %in% colnames(gene2go))) {
  stop("GENE2GO_FILE must contain columns: geneID, go_id")
}

gene2go <- gene2go %>%
  filter(!is.na(geneID), !is.na(go_id)) %>%
  distinct(geneID, go_id, .keep_all = TRUE)

universe_genes <- sort(unique(gene2go$geneID))
N <- length(universe_genes)

go_counts <- gene2go %>%
  distinct(geneID, go_id) %>%
  count(go_id, name = "go_universe_count")

go_meta <- gene2go %>%
  select(any_of(c("go_id", "go_name", "ontology"))) %>%
  distinct(go_id, .keep_all = TRUE)

run_hypergeom <- function(gene_set, set_name) {
  gene_set <- intersect(unique(gene_set), universe_genes)
  n <- length(gene_set)

  if (n == 0) {
    return(tibble())
  }

  gene2go %>%
    filter(geneID %in% gene_set) %>%
    distinct(geneID, go_id) %>%
    count(go_id, name = "overlap_count") %>%
    left_join(go_counts, by = "go_id") %>%
    mutate(
      set_name = set_name,
      set_size = n,
      universe_size = N,
      p_value = phyper(
        q = overlap_count - 1,
        m = go_universe_count,
        n = universe_size - go_universe_count,
        k = set_size,
        lower.tail = FALSE
      ),
      fold_enrichment = (overlap_count / set_size) / (go_universe_count / universe_size)
    ) %>%
    group_by(set_name) %>%
    mutate(
      p_adjust_bh = p.adjust(p_value, method = "BH"),
      p_adjust_bonferroni = p.adjust(p_value, method = "bonferroni")
    ) %>%
    ungroup() %>%
    left_join(go_meta, by = "go_id") %>%
    arrange(p_adjust_bonferroni, p_adjust_bh, p_value)
}

# Per-TF GO enrichment
by_tf <- targets %>%
  group_by(tf_geneID) %>%
  summarise(genes = list(unique(assigned_geneID)), .groups = "drop") %>%
  mutate(enrichment = map2(genes, tf_geneID, run_hypergeom)) %>%
  select(tf_geneID, enrichment) %>%
  unnest(enrichment)

if (nrow(by_tf) > 0) {
  write_tsv(by_tf, file.path(OUT_DIR, "go_enrichment_by_tf.tsv"))
}

# ERF clade-level target GO enrichment, if clade metadata is available
if (file.exists(ERF_CLADE_METADATA)) {
  clade <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE)
  if (all(c("geneID", "clade") %in% colnames(clade))) {
    erf_targets <- targets %>%
      inner_join(clade, by = c("tf_geneID" = "geneID"))

    by_clade <- erf_targets %>%
      group_by(clade) %>%
      summarise(genes = list(unique(assigned_geneID)), .groups = "drop") %>%
      mutate(set_name = paste0("ERF_clade_", clade)) %>%
      mutate(enrichment = map2(genes, set_name, run_hypergeom)) %>%
      select(clade, enrichment) %>%
      unnest(enrichment)

    if (nrow(by_clade) > 0) {
      write_tsv(by_clade, file.path(OUT_DIR, "go_enrichment_erf_clade_targets.tsv"))
    }
  }
}

message("GO enrichment complete.")
message("Done.")
