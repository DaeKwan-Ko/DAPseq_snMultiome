#!/usr/bin/env Rscript

# 01_make_txdb.R
# Build a sorghum TxDb object from the v3.1.1 GFF3 annotation.
# This is required because a prebuilt sorghum TxDb is not available.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(GenomicFeatures)
  library(txdbmaker)
})

if (!file.exists(GENOME_GFF3)) {
  stop("GENOME_GFF3 not found: ", GENOME_GFF3)
}

message("Building TxDb from GFF3: ", GENOME_GFF3)

txdb <- makeTxDbFromGFF(
  file = GENOME_GFF3,
  format = "gff3",
  dataSource = "Sorghum bicolor v3.1.1",
  organism = "Sorghum bicolor"
)

saveDb(txdb, file = TXDB_SQLITE)

message("Saved TxDb SQLite database to: ", TXDB_SQLITE)

# Save basic transcript/gene summary for QC
gene_ranges <- genes(txdb)
tx_summary <- data.frame(
  n_genes = length(gene_ranges),
  n_transcripts = length(transcripts(txdb)),
  n_exons = length(exons(txdb))
)

write.table(
  tx_summary,
  file = file.path(OUT_DIR, "txdb_summary.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

message("Done.")
