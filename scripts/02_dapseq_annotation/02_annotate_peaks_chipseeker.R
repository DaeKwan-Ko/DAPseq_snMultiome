#!/usr/bin/env Rscript

# 02_annotate_peaks_chipseeker.R
# Annotate DAP-seq peak BED files using ChIPseeker and the custom sorghum TxDb.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(ChIPseeker)
  library(GenomicFeatures)
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
})

if (!file.exists(TXDB_SQLITE)) {
  stop("TxDb SQLite file not found. Run 01_make_txdb.R first: ", TXDB_SQLITE)
}

if (!dir.exists(PEAK_BED_DIR)) {
  stop("PEAK_BED_DIR not found: ", PEAK_BED_DIR)
}

txdb <- loadDb(TXDB_SQLITE)

bed_files <- list.files(
  PEAK_BED_DIR,
  pattern = "_for_ChIPseeker\\.bed$|\\.bed$",
  full.names = TRUE
)

if (length(bed_files) == 0) {
  stop("No BED files found in PEAK_BED_DIR: ", PEAK_BED_DIR)
}

extract_gene_id_from_filename <- function(x) {
  # Expected filename starts with Sobic.xxxGxxxxxx
  id <- str_extract(basename(x), "Sobic\\.\\d{3}G\\d+")
  ifelse(is.na(id), tools::file_path_sans_ext(basename(x)), id)
}

annotate_one <- function(bed_file) {
  tf_id <- extract_gene_id_from_filename(bed_file)
  message("Annotating ", tf_id, ": ", basename(bed_file))

  peak <- readPeakFile(bed_file)

  anno <- annotatePeak(
    peak = peak,
    TxDb = txdb,
    tssRegion = TSS_REGION,
    verbose = FALSE
  )

  anno_df <- as.data.frame(anno)

  anno_df %>%
    mutate(
      tf_geneID = tf_id,
      source_file = basename(bed_file),
      peak_id = if ("name" %in% colnames(.)) name else paste(seqnames, start, end, sep = ":")
    ) %>%
    relocate(tf_geneID, source_file, peak_id)
}

all_annotations <- map_dfr(bed_files, annotate_one)

# Clean gene identifiers when possible.
# ChIPseeker usually stores the assigned gene in geneId.
if ("geneId" %in% colnames(all_annotations)) {
  all_annotations <- all_annotations %>%
    mutate(
      assigned_geneID = str_extract(as.character(geneId), "Sobic\\.\\d{3}G\\d+")
    )
} else {
  all_annotations$assigned_geneID <- NA_character_
}

# Keep a simplified annotation class.
all_annotations <- all_annotations %>%
  mutate(
    annotation_class = case_when(
      str_detect(annotation, regex("^Promoter", ignore_case = TRUE)) ~ "Promoter",
      str_detect(annotation, regex("^Exon", ignore_case = TRUE)) ~ "Exon",
      str_detect(annotation, regex("^Intron", ignore_case = TRUE)) ~ "Intron",
      str_detect(annotation, regex("^Downstream", ignore_case = TRUE)) ~ "Downstream",
      str_detect(annotation, regex("Distal Intergenic", ignore_case = TRUE)) ~ "Distal intergenic",
      TRUE ~ as.character(annotation)
    )
  )

out_file <- file.path(OUT_DIR, "all_peak_annotations.tsv")
write_tsv(all_annotations, out_file)

message("Saved full annotation table: ", out_file)
message("Number of annotated peaks: ", nrow(all_annotations))
message("Number of TFs annotated: ", dplyr::n_distinct(all_annotations$tf_geneID))
message("Done.")
