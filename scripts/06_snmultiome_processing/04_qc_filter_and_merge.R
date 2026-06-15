#!/usr/bin/env Rscript

# 04_qc_filter_and_merge.R
# Apply joint RNA/ATAC QC filters and merge samples.

source("scripts/06_snmultiome_processing/00_config.R")

suppressPackageStartupMessages({
  library(Seurat)
  library(Signac)
  library(dplyr)
  library(readr)
})

objects_file <- file.path(OUT_DIR, "snmultiome_sample_objects_raw.rds")
if (!file.exists(objects_file)) {
  stop("Raw sample object list not found. Run 03_load_cellranger_arc_seurat_signac.R first.")
}

objects <- readRDS(objects_file)

qc_one <- function(obj) {
  DefaultAssay(obj) <- "ATAC"

  # If TSS.enrichment and nucleosome_signal are absent, calculate where possible.
  if (!"nucleosome_signal" %in% colnames(obj@meta.data)) {
    obj <- NucleosomeSignal(obj)
  }

  if (!"TSS.enrichment" %in% colnames(obj@meta.data)) {
    obj <- TSSEnrichment(obj, fast = FALSE)
  }

  # Rename standard count columns if needed.
  if (!"nCount_RNA" %in% colnames(obj@meta.data)) stop("nCount_RNA missing")
  if (!"nCount_ATAC" %in% colnames(obj@meta.data)) stop("nCount_ATAC missing")

  obj$pass_joint_qc <- with(
    obj@meta.data,
    nCount_RNA >= MIN_RNA_COUNTS &
      nCount_RNA <= MAX_RNA_COUNTS &
      nCount_ATAC >= MIN_ATAC_COUNTS &
      nCount_ATAC <= MAX_ATAC_COUNTS &
      TSS.enrichment > MIN_TSS_ENRICHMENT &
      nucleosome_signal < MAX_NUCLEOSOME_SIGNAL
  )

  obj <- subset(obj, subset = pass_joint_qc)
  obj
}

filtered <- lapply(objects, qc_one)

qc_summary <- bind_rows(lapply(names(objects), function(sample_id) {
  before <- ncol(objects[[sample_id]])
  after <- ncol(filtered[[sample_id]])
  data.frame(
    sample_id = sample_id,
    nuclei_before_qc = before,
    nuclei_after_qc = after,
    nuclei_removed = before - after
  )
}))

write_tsv(qc_summary, file.path(OUT_DIR, "joint_qc_summary_by_sample.tsv"))

merged <- merge(
  filtered[[1]],
  y = filtered[-1],
  add.cell.ids = names(filtered),
  project = "sorghum_snmultiome"
)

saveRDS(merged, MERGED_OBJECT)

message("Saved merged QC-filtered object: ", MERGED_OBJECT)
message("Nuclei after QC: ", ncol(merged))
print(qc_summary)
