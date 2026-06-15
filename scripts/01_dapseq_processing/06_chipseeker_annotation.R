#!/usr/bin/env Rscript
# Annotate DAP-seq peaks with ChIPseeker using a custom sorghum TxDb.
# Usage:
#   Rscript 06_chipseeker_annotation.R \
#     --bed_dir /path/to/peaks_bed_for_chipseeker \
#     --gff3 /path/to/Sbicolor_454_v3.1.1.gene.gff3 \
#     --out_dir /path/to/chipseeker_annotation

suppressPackageStartupMessages({
  library(optparse)
  library(ChIPseeker)
  library(GenomicFeatures)
  library(rtracklayer)
  library(dplyr)
  library(readr)
  library(stringr)
})

option_list <- list(
  make_option("--bed_dir", type = "character", help = "Directory containing BED files."),
  make_option("--gff3", type = "character", help = "Sorghum GFF3 annotation file."),
  make_option("--out_dir", type = "character", help = "Output directory."),
  make_option("--promoter_upstream", type = "integer", default = 5000, help = "Promoter upstream distance [default %default]."),
  make_option("--promoter_downstream", type = "integer", default = 0, help = "Promoter downstream distance [default %default].")
)
opt <- parse_args(OptionParser(option_list = option_list))

if (is.null(opt$bed_dir) || is.null(opt$gff3) || is.null(opt$out_dir)) {
  stop("Required arguments: --bed_dir, --gff3, --out_dir")
}

if (!dir.exists(opt$bed_dir)) stop("BED directory does not exist: ", opt$bed_dir)
if (!file.exists(opt$gff3)) stop("GFF3 file does not exist: ", opt$gff3)
dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "annotated_peaks"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(opt$out_dir, "candidate_targets"), recursive = TRUE, showWarnings = FALSE)

message("Building TxDb from GFF3: ", opt$gff3)
txdb <- makeTxDbFromGFF(opt$gff3, format = "gff3")

bed_files <- list.files(opt$bed_dir, pattern = "\\.bed$", full.names = TRUE)
if (length(bed_files) == 0) stop("No BED files found in: ", opt$bed_dir)

is_candidate_target_annotation <- function(annotation) {
  ann <- as.character(annotation)
  str_detect(ann, "^Promoter") |
    str_detect(ann, "Exon") |
    str_detect(ann, "Intron") |
    str_detect(ann, "Downstream \\(<=300")
}

summary_list <- vector("list", length(bed_files))

for (i in seq_along(bed_files)) {
  bed <- bed_files[[i]]
  gene_id <- tools::file_path_sans_ext(basename(bed))
  message("Annotating ", gene_id, " (", i, "/", length(bed_files), ")")

  peak_gr <- readPeakFile(bed)
  anno <- annotatePeak(
    peak_gr,
    TxDb = txdb,
    tssRegion = c(-opt$promoter_upstream, opt$promoter_downstream),
    verbose = FALSE
  )

  anno_df <- as.data.frame(anno)
  anno_df <- anno_df %>% mutate(TF_geneID = gene_id, .before = 1)

  annotated_out <- file.path(opt$out_dir, "annotated_peaks", paste0(gene_id, ".annotated_peaks.tsv"))
  write_tsv(anno_df, annotated_out)

  target_df <- anno_df %>%
    filter(is_candidate_target_annotation(annotation)) %>%
    filter(!is.na(geneId)) %>%
    distinct(TF_geneID, geneId, annotation, distanceToTSS)

  target_out <- file.path(opt$out_dir, "candidate_targets", paste0(gene_id, ".candidate_targets.tsv"))
  write_tsv(target_df, target_out)

  summary_list[[i]] <- tibble(
    TF_geneID = gene_id,
    total_peaks = nrow(anno_df),
    promoter_peaks = sum(str_detect(as.character(anno_df$annotation), "^Promoter"), na.rm = TRUE),
    genic_peaks = sum(str_detect(as.character(anno_df$annotation), "Exon|Intron"), na.rm = TRUE),
    downstream_300bp_peaks = sum(str_detect(as.character(anno_df$annotation), "Downstream \\(<=300"), na.rm = TRUE),
    distal_intergenic_peaks = sum(str_detect(as.character(anno_df$annotation), "Distal Intergenic"), na.rm = TRUE),
    candidate_target_gene_count = n_distinct(target_df$geneId)
  )
}

summary_df <- bind_rows(summary_list)
write_tsv(summary_df, file.path(opt$out_dir, "chipseeker_annotation_summary.tsv"))

all_targets <- list.files(file.path(opt$out_dir, "candidate_targets"), pattern = "candidate_targets.tsv$", full.names = TRUE) %>%
  lapply(read_tsv, show_col_types = FALSE) %>%
  bind_rows()
write_tsv(all_targets, file.path(opt$out_dir, "all_TF_candidate_targets.tsv"))

message("Done. Results written to: ", opt$out_dir)
