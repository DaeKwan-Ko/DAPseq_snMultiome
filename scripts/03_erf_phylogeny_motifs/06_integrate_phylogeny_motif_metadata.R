#!/usr/bin/env Rscript

# 06_integrate_phylogeny_motif_metadata.R
# Combine ERF metadata, peak counts, beta2 signatures and top motif summaries.

source("scripts/03_erf_phylogeny_motifs/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

meta_file <- file.path(OUT_DIR, "erf_metadata_prepared.tsv")
if (!file.exists(meta_file)) {
  stop("ERF metadata not found. Run 01_prepare_erf_metadata.R first: ", meta_file)
}

meta <- read_tsv(meta_file, show_col_types = FALSE)

# Optional beta2 signature table
beta2_file <- file.path(OUT_DIR, "erf_beta2_signature.tsv")
if (file.exists(beta2_file)) {
  beta2 <- read_tsv(beta2_file, show_col_types = FALSE)
  meta <- meta %>% left_join(beta2, by = "geneID")
} else {
  meta <- meta %>% mutate(beta2_signature = NA_character_, beta2_signature_class = NA_character_)
}

# Optional MEME top motif table
motif_file <- file.path(OUT_DIR, "meme_top_motif_summary.tsv")
if (file.exists(motif_file)) {
  motifs <- read_tsv(motif_file, show_col_types = FALSE)
  meta <- meta %>%
    left_join(
      motifs %>% select(geneID, motif_id, motif_rank, width, nsites, e_value, top_motif_label),
      by = "geneID"
    )
} else {
  meta <- meta %>% mutate(motif_id = NA_character_, e_value = NA_real_)
}

# Optional peak-count table. Try common column names.
if (file.exists(PEAK_SUMMARY)) {
  peaks <- read_tsv(PEAK_SUMMARY, show_col_types = FALSE)

  if ("tf_geneID" %in% colnames(peaks)) {
    peaks <- peaks %>% rename(geneID = tf_geneID)
  }
  if ("peak_count_chr1_chr10" %in% colnames(peaks) && !("n_peaks" %in% colnames(peaks))) {
    peaks <- peaks %>% rename(n_peaks = peak_count_chr1_chr10)
  }

  keep_cols <- intersect(c("geneID", "n_peaks", "bound_gene_count", "genome_coverage_bp"), colnames(peaks))
  if (all(c("geneID") %in% keep_cols)) {
    meta <- meta %>% left_join(peaks %>% select(all_of(keep_cols)), by = "geneID")
  }
}

meta <- meta %>%
  arrange(clade, geneID)

out_file <- file.path(OUT_DIR, "erf_phylogeny_motif_annotation_table.tsv")
write_tsv(meta, out_file)

message("Saved integrated ERF annotation table: ", out_file)
message("Done.")
