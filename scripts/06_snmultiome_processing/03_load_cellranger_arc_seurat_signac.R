#!/usr/bin/env Rscript

# 03_load_cellranger_arc_seurat_signac.R
# Load Cell Ranger ARC outputs into Seurat/Signac objects.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(GenomicRanges)
  library(GenomeInfoDb)
  library(readr)
  library(dplyr)
  library(stringr)
  library(future)
})

plan("sequential")

if (!file.exists(SAMPLE_METADATA)) {
  stop("Sample metadata not found: ", SAMPLE_METADATA)
}

sample_meta <- read_tsv(SAMPLE_METADATA, show_col_types = FALSE)
if (!"sample_id" %in% colnames(sample_meta)) {
  stop("sample_metadata_snmultiome.tsv must contain a sample_id column")
}

load_one_arc_sample <- function(sample_row) {
  sample_id <- sample_row$sample_id
  outs_dir <- file.path(CELLRANGER_ARC_DIR, sample_id, "outs")

  h5_file <- file.path(outs_dir, "filtered_feature_bc_matrix.h5")
  fragments_file <- file.path(outs_dir, "atac_fragments.tsv.gz")
  per_barcode_metrics <- file.path(outs_dir, "per_barcode_metrics.csv")

  if (!file.exists(h5_file)) stop("Missing filtered_feature_bc_matrix.h5 for ", sample_id)
  if (!file.exists(fragments_file)) stop("Missing atac_fragments.tsv.gz for ", sample_id)

  message("Loading sample: ", sample_id)
  counts <- Read10X_h5(h5_file)

  # Cell Ranger ARC usually returns named assays: Gene Expression and Peaks.
  rna_counts <- counts[["Gene Expression"]]
  atac_counts <- counts[["Peaks"]]

  obj <- CreateSeuratObject(
    counts = rna_counts,
    assay = "RNA",
    project = sample_id,
    meta.data = data.frame(sample_id = sample_id, row.names = colnames(rna_counts))
  )

  chrom_assay <- CreateChromatinAssay(
    counts = atac_counts,
    sep = c(":", "-"),
    fragments = fragments_file,
    min.cells = 1,
    min.features = 1
  )

  obj[["ATAC"]] <- chrom_assay
  obj$condition <- sample_row$condition
  obj$treatment <- sample_row$treatment
  obj$genotype <- sample_row$genotype
  obj$species <- sample_row$species
  obj$tissue <- sample_row$tissue

  # Add Cell Ranger per-barcode metrics if present.
  if (file.exists(per_barcode_metrics)) {
    metrics <- read_csv(per_barcode_metrics, show_col_types = FALSE)
    barcode_col <- intersect(c("barcode", "Barcode"), colnames(metrics))[1]
    if (!is.na(barcode_col)) {
      metrics <- metrics %>% as.data.frame()
      rownames(metrics) <- metrics[[barcode_col]]
      common <- intersect(colnames(obj), rownames(metrics))
      obj <- AddMetaData(obj, metadata = metrics[common, , drop = FALSE])
    }
  }

  obj
}

objects <- lapply(seq_len(nrow(sample_meta)), function(i) load_one_arc_sample(sample_meta[i, ]))
names(objects) <- sample_meta$sample_id

out_file <- file.path(OUT_DIR, "snmultiome_sample_objects_raw.rds")
saveRDS(objects, out_file)

message("Saved raw sample object list: ", out_file)
message("Loaded samples: ", paste(names(objects), collapse = ", "))
