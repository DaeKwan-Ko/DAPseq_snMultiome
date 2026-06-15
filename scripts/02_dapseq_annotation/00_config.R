# 00_config.R
# Central configuration file for DAP-seq peak annotation.
# Edit paths here before running the scripts in this folder.

suppressPackageStartupMessages({
  library(fs)
})

# -------------------------------------------------------------------------
# Project paths
# -------------------------------------------------------------------------

# Root of the GitHub repository after cloning.
# You can also run from the repository root and keep this default.
PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)

# Input metadata
TF_METADATA <- file.path(PROJECT_DIR, "metadata", "tf_list_142.tsv")
ERF_CLADE_METADATA <- file.path(PROJECT_DIR, "metadata", "erf_clade_annotation.tsv")

# Optional GO annotation table.
# Required columns: geneID, go_id
# Optional columns: go_name, ontology
GENE2GO_FILE <- file.path(PROJECT_DIR, "metadata", "sorghum_gene2go.tsv")

# Reference files
GENOME_GFF3 <- file.path(PROJECT_DIR, "reference", "Sbicolor_454_v3.1.1.gene.gff3")

# Input BED directory from scripts/01_dapseq_processing.
# Expected files: *_for_ChIPseeker.bed
PEAK_BED_DIR <- file.path(PROJECT_DIR, "results", "01_dapseq_processing", "peaks_for_chipseeker")

# Output directory
OUT_DIR <- file.path(PROJECT_DIR, "results", "02_dapseq_annotation")
TXDB_SQLITE <- file.path(OUT_DIR, "sorghum_v3.1.1.txdb.sqlite")

# -------------------------------------------------------------------------
# Annotation parameters
# -------------------------------------------------------------------------

PROMOTER_UPSTREAM_BP <- 5000
PROMOTER_DOWNSTREAM_BP <- 0
DOWNSTREAM_MAX_BP <- 300

# ChIPseeker will annotate peaks using this promoter window.
TSS_REGION <- c(-PROMOTER_UPSTREAM_BP, PROMOTER_DOWNSTREAM_BP)

# A peak-associated gene is retained as a conservative candidate target if
# the ChIPseeker annotation is promoter-proximal, genic, or <=300 bp downstream.
CONSERVATIVE_ANNOTATION_CLASSES <- c(
  "Promoter",
  "Exon",
  "Intron",
  "Downstream"
)

# -------------------------------------------------------------------------
# Utility functions
# -------------------------------------------------------------------------

dir_create(OUT_DIR, recurse = TRUE)
dir_create(file.path(OUT_DIR, "plots"), recurse = TRUE)
dir_create(file.path(OUT_DIR, "intermediate"), recurse = TRUE)

message("PROJECT_DIR: ", PROJECT_DIR)
message("OUT_DIR: ", OUT_DIR)
