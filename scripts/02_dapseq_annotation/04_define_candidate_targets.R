#!/usr/bin/env Rscript

# 04_define_candidate_targets.R
# Define conservative DAP-seq-associated candidate target genes.
#
# Candidate target definition:
# - promoter-proximal peaks within 5 kb upstream of TSS;
# - genic peaks;
# - downstream peaks within 300 bp of annotated genes.
# Distal intergenic peaks are excluded.

source("scripts/02_dapseq_annotation/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidyr)
})

annotation_file <- file.path(OUT_DIR, "all_peak_annotations.tsv")
if (!file.exists(annotation_file)) {
  stop("Annotation file not found. Run 02_annotate_peaks_chipseeker.R first: ", annotation_file)
}

anno <- read_tsv(annotation_file, show_col_types = FALSE)

required <- c("tf_geneID", "assigned_geneID", "annotation_class")
missing <- setdiff(required, colnames(anno))
if (length(missing) > 0) {
  stop("Missing required columns in annotation table: ", paste(missing, collapse = ", "))
}

# ChIPseeker distanceToTSS is useful for downstream filtering.
# Downstream peaks can have positive or negative values depending on strand/context,
# so here we keep ChIPseeker Downstream class and require absolute distance <= 300
# when distanceToTSS is available. If unavailable, all Downstream-class peaks are kept.
candidate_peak_gene <- anno %>%
  filter(!is.na(assigned_geneID)) %>%
  mutate(
    keep_candidate = case_when(
      annotation_class %in% c("Promoter", "Exon", "Intron") ~ TRUE,
      annotation_class == "Downstream" & "distanceToTSS" %in% colnames(.) ~ abs(distanceToTSS) <= DOWNSTREAM_MAX_BP,
      annotation_class == "Downstream" ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  filter(keep_candidate) %>%
  select(
    tf_geneID,
    assigned_geneID,
    peak_id,
    seqnames,
    start,
    end,
    annotation,
    annotation_class,
    distanceToTSS,
    everything()
  )

candidate_targets_long <- candidate_peak_gene %>%
  distinct(tf_geneID, assigned_geneID, annotation_class) %>%
  arrange(tf_geneID, assigned_geneID)

candidate_targets_by_tf <- candidate_targets_long %>%
  group_by(tf_geneID) %>%
  summarise(
    n_candidate_targets = n_distinct(assigned_geneID),
    candidate_targets = paste(sort(unique(assigned_geneID)), collapse = ";"),
    .groups = "drop"
  )

target_count_by_tf <- candidate_peak_gene %>%
  group_by(tf_geneID, assigned_geneID) %>%
  summarise(
    n_supporting_peaks = n_distinct(peak_id),
    annotation_classes = paste(sort(unique(annotation_class)), collapse = ";"),
    min_abs_distance_to_tss = if ("distanceToTSS" %in% colnames(.)) min(abs(distanceToTSS), na.rm = TRUE) else NA_real_,
    .groups = "drop"
  ) %>%
  arrange(tf_geneID, assigned_geneID)

write_tsv(candidate_peak_gene, file.path(OUT_DIR, "dapseq_candidate_peak_gene_links.tsv"))
write_tsv(candidate_targets_long, file.path(OUT_DIR, "dapseq_candidate_targets_long.tsv"))
write_tsv(candidate_targets_by_tf, file.path(OUT_DIR, "dapseq_candidate_targets_by_tf.tsv"))
write_tsv(target_count_by_tf, file.path(OUT_DIR, "dapseq_target_count_by_tf.tsv"))

message("Saved conservative candidate target tables.")
message("Number of TFs with candidate targets: ", n_distinct(candidate_targets_long$tf_geneID))
message("Number of unique candidate target genes: ", n_distinct(candidate_targets_long$assigned_geneID))
message("Done.")
