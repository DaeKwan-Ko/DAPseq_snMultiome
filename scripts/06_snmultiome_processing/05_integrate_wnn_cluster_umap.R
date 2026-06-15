#!/usr/bin/env Rscript

# 05_integrate_wnn_cluster_umap.R
# Normalize RNA and ATAC assays and perform WNN integration, clustering and UMAP.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
})

if (!file.exists(MERGED_OBJECT)) {
  stop("Merged object not found. Run 04_qc_filter_and_merge.R first: ", MERGED_OBJECT)
}

obj <- readRDS(MERGED_OBJECT)

# RNA processing
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE)
obj <- RunPCA(obj, assay = "SCT", verbose = FALSE)

# ATAC processing
DefaultAssay(obj) <- "ATAC"
obj <- FindTopFeatures(obj, min.cutoff = "q0")
obj <- RunTFIDF(obj)
obj <- RunSVD(obj)

# Weighted nearest-neighbor integration
obj <- FindMultiModalNeighbors(
  obj,
  reduction.list = list("pca", "lsi"),
  dims.list = list(RNA_DIMS, ATAC_DIMS)
)

obj <- RunUMAP(
  obj,
  nn.name = "weighted.nn",
  assay = "RNA",
  reduction.name = "wnn.umap",
  reduction.key = "wnnUMAP_"
)

obj <- FindClusters(
  obj,
  graph.name = "wsnn",
  algorithm = 3,
  resolution = CLUSTER_RESOLUTION,
  verbose = FALSE
)

saveRDS(obj, WNN_OBJECT)

cluster_summary <- obj@meta.data %>%
  count(seurat_clusters, condition, name = "n_nuclei") %>%
  arrange(seurat_clusters, condition)

write_tsv(cluster_summary, file.path(OUT_DIR, "cluster_condition_nucleus_counts.tsv"))

message("Saved WNN-integrated object: ", WNN_OBJECT)
message("Clusters:")
print(table(obj$seurat_clusters))
