#!/usr/bin/env Rscript

# 07_calculate_gene_activity.R
# Calculate ATAC-derived gene activity matrix using Signac.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicFeatures)
  library(GenomicRanges)
  library(dplyr)
  library(readr)
})

if (!file.exists(ANNOTATED_OBJECT)) {
  stop("Annotated object not found. Run 06_marker_based_celltype_annotation.R first: ", ANNOTATED_OBJECT)
}

obj <- readRDS(ANNOTATED_OBJECT)

if (!file.exists(GENOME_GFF3)) {
  message("GENOME_GFF3 not found: ", GENOME_GFF3)
  message("Skipping GeneActivity. Provide reference annotation and rerun.")
  quit(save = "no", status = 0)
}

message("Building gene annotation from GFF3: ", GENOME_GFF3)
txdb <- makeTxDbFromGFF(GENOME_GFF3, format = "gff3")
gene_ranges <- genes(txdb)

# Assign annotation to ATAC assay.
DefaultAssay(obj) <- "ATAC"
Annotation(obj) <- gene_ranges

message("Calculating gene activity matrix...")
gene.activities <- GeneActivity(obj)

obj[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(obj) <- "ACTIVITY"
obj <- NormalizeData(obj, assay = "ACTIVITY")
obj <- ScaleData(obj, assay = "ACTIVITY", features = rownames(obj[["ACTIVITY"]]), verbose = FALSE)

saveRDS(obj, GENE_ACTIVITY_OBJECT)
saveRDS(gene.activities, file.path(OUT_DIR, "gene_activity_matrix.rds"))

message("Saved object with gene activity assay: ", GENE_ACTIVITY_OBJECT)
message("Gene activity features: ", nrow(gene.activities))
