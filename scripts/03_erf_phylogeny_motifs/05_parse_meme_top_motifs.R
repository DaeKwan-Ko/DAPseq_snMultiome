#!/usr/bin/env Rscript

# 05_parse_meme_top_motifs.R
# Parse MEME outputs and extract top motif summaries for each ERF.
#
# This script recursively searches MEME_DIR for meme.txt and .meme files.
# It extracts motif ID, width, site count and E-value when available.

source("scripts/03_erf_phylogeny_motifs/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

if (!dir.exists(MEME_DIR)) {
  message("MEME_DIR not found: ", MEME_DIR)
  message("Skipping motif parsing.")
  quit(save = "no", status = 0)
}

meme_files <- list.files(
  MEME_DIR,
  pattern = "(meme\\.txt$|\\.meme$)",
  recursive = TRUE,
  full.names = TRUE
)

if (length(meme_files) == 0) {
  message("No MEME files found in: ", MEME_DIR)
  quit(save = "no", status = 0)
}

extract_gene_id <- function(path) {
  id <- str_extract(path, "Sobic\\.\\d{3}G\\d+")
  ifelse(is.na(id), basename(dirname(path)), id)
}

parse_one_meme <- function(path) {
  lines <- readLines(path, warn = FALSE)
  motif_lines <- grep("^MOTIF\\s+", lines, value = TRUE)

  if (length(motif_lines) == 0) {
    return(tibble())
  }

  map_dfr(seq_along(motif_lines), function(i) {
    line <- motif_lines[[i]]
    motif_id <- str_split(line, "\\s+")[[1]][2]

    # Find metadata line following the motif line.
    idx <- grep(paste0("^", fixed(line)), lines, fixed = TRUE)[1]
    window <- lines[idx:min(length(lines), idx + 10)]
    meta <- window[str_detect(window, "letter-probability matrix")]

    width <- NA_integer_
    nsites <- NA_integer_
    evalue <- NA_real_

    if (length(meta) > 0) {
      width_chr <- str_match(meta[1], "w=\\s*([0-9]+)")[,2]
      nsites_chr <- str_match(meta[1], "nsites=\\s*([0-9]+)")[,2]
      evalue_chr <- str_match(meta[1], "E=\\s*([0-9.eE+-]+)")[,2]
      width <- suppressWarnings(as.integer(width_chr))
      nsites <- suppressWarnings(as.integer(nsites_chr))
      evalue <- suppressWarnings(as.numeric(evalue_chr))
    }

    tibble(
      geneID = extract_gene_id(path),
      meme_file = path,
      motif_rank = i,
      motif_id = motif_id,
      width = width,
      nsites = nsites,
      e_value = evalue
    )
  })
}

motifs <- map_dfr(meme_files, parse_one_meme)

if (nrow(motifs) == 0) {
  message("No motif records parsed.")
  quit(save = "no", status = 0)
}

top_motifs <- motifs %>%
  group_by(geneID) %>%
  arrange(e_value, motif_rank, .by_group = TRUE) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    top_motif_label = paste0(motif_id, "_w", width)
  )

write_tsv(motifs, file.path(OUT_DIR, "meme_all_motif_summary.tsv"))
write_tsv(top_motifs, file.path(OUT_DIR, "meme_top_motif_summary.tsv"))

message("Saved MEME motif summaries.")
message("Number of MEME files parsed: ", length(meme_files))
message("Number of genes with top motif: ", n_distinct(top_motifs$geneID))
message("Done.")
