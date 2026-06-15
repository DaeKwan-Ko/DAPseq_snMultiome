#!/usr/bin/env Rscript

# 08_qc_check_erf_outputs.R
# Quick QC check for ERF phylogeny and motif outputs.

source("scripts/03_erf_phylogeny_motifs/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(Biostrings)
})

check_file <- function(path) {
  data.frame(
    file = basename(path),
    path = path,
    exists = file.exists(path),
    size_bytes = ifelse(file.exists(path), file.info(path)$size, NA_real_)
  )
}

paths <- c(
  file.path(OUT_DIR, "erf_metadata_prepared.tsv"),
  file.path(OUT_DIR, "sorghum_erf_ap2_domains.fa"),
  file.path(OUT_DIR, "sorghum_erf_ap2_domains_aligned.fa"),
  file.path(OUT_DIR, "sorghum_erf_dbd_tree.nwk"),
  file.path(OUT_DIR, "erf_beta2_signature.tsv"),
  file.path(OUT_DIR, "meme_top_motif_summary.tsv"),
  file.path(OUT_DIR, "erf_phylogeny_motif_annotation_table.tsv"),
  file.path(OUT_DIR, "ceqlogo_commands.sh")
)

qc <- bind_rows(lapply(paths, check_file))
write_tsv(qc, file.path(OUT_DIR, "erf_phylogeny_motif_output_qc.tsv"))
print(qc)

if (file.exists(file.path(OUT_DIR, "erf_metadata_prepared.tsv"))) {
  meta <- read_tsv(file.path(OUT_DIR, "erf_metadata_prepared.tsv"), show_col_types = FALSE)
  message("ERF metadata rows: ", nrow(meta))
  if ("clade" %in% colnames(meta)) {
    print(table(meta$clade, useNA = "ifany"))
  }
}

if (file.exists(file.path(OUT_DIR, "sorghum_erf_ap2_domains.fa"))) {
  seqs <- readAAStringSet(file.path(OUT_DIR, "sorghum_erf_ap2_domains.fa"))
  message("Extracted DBD sequences: ", length(seqs))
}

message("Saved QC table: ", file.path(OUT_DIR, "erf_phylogeny_motif_output_qc.tsv"))
