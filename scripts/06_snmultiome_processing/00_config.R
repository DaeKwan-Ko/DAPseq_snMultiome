# 00_config.R
# Central R configuration for sorghum root snMultiome processing.

suppressPackageStartupMessages({
  library(fs)
})

PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)

METADATA_DIR <- file.path(PROJECT_DIR, "metadata")
OUT_DIR <- file.path(PROJECT_DIR, "results", "06_snmultiome_processing")
INPUT_DIR <- file.path(OUT_DIR, "input")
CELLRANGER_ARC_DIR <- file.path(INPUT_DIR, "cellranger_arc")

SAMPLE_METADATA <- file.path(METADATA_DIR, "sample_metadata_snmultiome.tsv")
TF_METADATA <- file.path(METADATA_DIR, "tf_list_142.tsv")
ERF_CLADE_METADATA <- file.path(METADATA_DIR, "erf_clade_annotation.tsv")
SALINITY_ERF_METADATA <- file.path(METADATA_DIR, "salinity_linked_erfs.tsv")
ROOT_MARKER_GENES <- file.path(METADATA_DIR, "root_marker_genes.tsv")
CANDIDATE_TARGET_GENES <- file.path(METADATA_DIR, "candidate_target_genes.tsv")

# Reference annotation used for gene activity.
GENOME_GFF3 <- file.path(PROJECT_DIR, "reference", "Sbicolor_454_v3.1.1.gene.gff3")

# Output objects
MERGED_OBJECT <- file.path(OUT_DIR, "snmultiome_merged_qc_filtered.rds")
WNN_OBJECT <- file.path(OUT_DIR, "snmultiome_wnn_integrated.rds")
ANNOTATED_OBJECT <- file.path(OUT_DIR, "snmultiome_annotated.rds")
GENE_ACTIVITY_OBJECT <- file.path(OUT_DIR, "snmultiome_with_gene_activity.rds")

# QC thresholds used in the manuscript workflow.
MIN_RNA_COUNTS <- 200
MAX_RNA_COUNTS <- 50000
MIN_ATAC_COUNTS <- 200
MAX_ATAC_COUNTS <- 50000
MIN_TSS_ENRICHMENT <- 1.2
MAX_NUCLEOSOME_SIGNAL <- 0.35

# Dimensionality settings
RNA_DIMS <- 1:30
ATAC_DIMS <- 2:30
WNN_DIMS <- 1:30
CLUSTER_RESOLUTION <- 0.5

dir_create(OUT_DIR, recurse = TRUE)
dir_create(INPUT_DIR, recurse = TRUE)
dir_create(CELLRANGER_ARC_DIR, recurse = TRUE)
dir_create(file.path(OUT_DIR, "plots"), recurse = TRUE)
dir_create(file.path(OUT_DIR, "tables"), recurse = TRUE)

message("PROJECT_DIR: ", PROJECT_DIR)
message("OUT_DIR: ", OUT_DIR)
