#!/usr/bin/env Rscript

# 06_marker_based_celltype_annotation.R
# Annotate root cell populations using marker gene modules.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
})

if (!file.exists(WNN_OBJECT)) {
  stop("WNN object not found. Run 05_integrate_wnn_cluster_umap.R first: ", WNN_OBJECT)
}

obj <- readRDS(WNN_OBJECT)

if (!file.exists(ROOT_MARKER_GENES)) {
  message("Marker gene table not found: ", ROOT_MARKER_GENES)
  message("Creating placeholder annotation by cluster. Replace with marker-based labels after preparing metadata/root_marker_genes.tsv.")
  obj$celltype <- paste0("Cluster_", obj$seurat_clusters)
  saveRDS(obj, ANNOTATED_OBJECT)
  write_tsv(
    obj@meta.data %>%
      tibble::rownames_to_column("barcode") %>%
      select(barcode, sample_id, condition, seurat_clusters, celltype),
    file.path(OUT_DIR, "celltype_annotation_metadata.tsv")
  )
  quit(save = "no", status = 0)
}

markers <- read_tsv(ROOT_MARKER_GENES, show_col_types = FALSE)
if (!all(c("celltype", "geneID") %in% colnames(markers))) {
  stop("root_marker_genes.tsv must contain columns: celltype, geneID")
}

DefaultAssay(obj) <- "SCT"

marker_list <- split(markers$geneID, markers$celltype)
marker_list <- lapply(marker_list, function(x) intersect(unique(x), rownames(obj)))

# Remove empty marker sets.
marker_list <- marker_list[lengths(marker_list) > 0]

if (length(marker_list) == 0) {
  stop("No marker genes were found in the object. Check gene identifiers.")
}

obj <- AddModuleScore(
  object = obj,
  features = marker_list,
  name = "markerModule_",
  assay = "SCT"
)

score_cols <- paste0("markerModule_", seq_along(marker_list))
names(score_cols) <- names(marker_list)

# Average module scores by cluster.
score_df <- obj@meta.data %>%
  select(seurat_clusters, all_of(score_cols)) %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean, na.rm = TRUE), .groups = "drop")

score_long <- score_df %>%
  pivot_longer(
    cols = all_of(score_cols),
    names_to = "score_col",
    values_to = "mean_module_score"
  ) %>%
  mutate(celltype_candidate = names(score_cols)[match(score_col, score_cols)])

cluster_labels <- score_long %>%
  group_by(seurat_clusters) %>%
  slice_max(mean_module_score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  transmute(seurat_clusters, celltype = celltype_candidate)

obj@meta.data <- obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  left_join(cluster_labels, by = "seurat_clusters") %>%
  tibble::column_to_rownames("barcode")

obj$celltype <- ifelse(is.na(obj$celltype), "unassigned", obj$celltype)

saveRDS(obj, ANNOTATED_OBJECT)

write_tsv(score_long, file.path(OUT_DIR, "marker_module_scores_by_cluster.tsv"))
write_tsv(cluster_labels, file.path(OUT_DIR, "cluster_to_celltype_annotation.tsv"))

cell_meta <- obj@meta.data %>%
  tibble::rownames_to_column("barcode") %>%
  select(barcode, sample_id, condition, seurat_clusters, celltype, nCount_RNA, nCount_ATAC, TSS.enrichment, nucleosome_signal)

write_tsv(cell_meta, file.path(OUT_DIR, "celltype_annotation_metadata.tsv"))

counts <- obj@meta.data %>%
  count(celltype, condition, name = "n_nuclei") %>%
  group_by(condition) %>%
  mutate(fraction_of_condition = n_nuclei / sum(n_nuclei)) %>%
  ungroup()

write_tsv(counts, file.path(OUT_DIR, "celltype_condition_nucleus_counts.tsv"))

message("Saved annotated object: ", ANNOTATED_OBJECT)
print(counts)
