# 00_config.R
# Central configuration for key plotting scripts.

suppressPackageStartupMessages({
  library(fs)
})

PROJECT_DIR <- normalizePath(getwd(), mustWork = FALSE)

METADATA_DIR <- file.path(PROJECT_DIR, "metadata")
OUT_DIR <- file.path(PROJECT_DIR, "results", "08_key_plotting")
PLOT_DIR <- file.path(OUT_DIR, "plots")
SOURCE_DATA_DIR <- file.path(OUT_DIR, "source_data")

# Metadata
TF_METADATA <- file.path(METADATA_DIR, "tf_list_142.tsv")
ERF_CLADE_METADATA <- file.path(METADATA_DIR, "erf_clade_annotation.tsv")
SALINITY_ERF_METADATA <- file.path(METADATA_DIR, "salinity_linked_erfs.tsv")
SPECIES_METADATA <- file.path(METADATA_DIR, "species_list_beta2_signature.tsv")

# Upstream result directories
DAPSEQ_ANNOT_DIR <- file.path(PROJECT_DIR, "results", "02_dapseq_annotation")
ERF_PHYLO_DIR <- file.path(PROJECT_DIR, "results", "03_erf_phylogeny_motifs")
EVOLUTION_DIR <- file.path(PROJECT_DIR, "results", "04_evolutionary_analysis")
SNMULTIOME_DIR <- file.path(PROJECT_DIR, "results", "06_snmultiome_processing")
NETWORK_DIR <- file.path(PROJECT_DIR, "results", "07_peak_to_gene_networks")

# Upstream files
PEAK_ANNOT_SUMMARY <- file.path(DAPSEQ_ANNOT_DIR, "peak_annotation_summary_by_tf.tsv")
ERF_TARGET_CLADE_COUNTS <- file.path(DAPSEQ_ANNOT_DIR, "erf_target_clade_counts.tsv")
ERF_CLADE_BIASED_TARGETS <- file.path(DAPSEQ_ANNOT_DIR, "erf_clade_biased_recurrent_targets.tsv")

ERF_PHYLO_MOTIF_TABLE <- file.path(ERF_PHYLO_DIR, "erf_phylogeny_motif_annotation_table.tsv")
ERF_TREE_NWK <- file.path(ERF_PHYLO_DIR, "sorghum_erf_dbd_tree.nwk")
MEME_TOP_MOTIFS <- file.path(ERF_PHYLO_DIR, "meme_top_motif_summary.tsv")

KAKS_PAIRS <- file.path(EVOLUTION_DIR, "kaks_ortholog_pairs_annotated.tsv")
KAKS_SUMMARY <- file.path(EVOLUTION_DIR, "kaks_summary_by_clade.tsv")
KAKS_TESTS <- file.path(EVOLUTION_DIR, "kaks_wilcoxon_tests.tsv")
BETA2_LINEAGE_SUMMARY <- file.path(EVOLUTION_DIR, "cross_lineage_beta2_summary_for_plotting.tsv")

CELLTYPE_COUNTS <- file.path(SNMULTIOME_DIR, "celltype_condition_nucleus_counts.tsv")
ERF_RNA_ACTIVITY_FC <- file.path(SNMULTIOME_DIR, "erf_rna_activity_log2fc_by_celltype.tsv")
TARGET_RNA_ACTIVITY_FC <- file.path(SNMULTIOME_DIR, "target_gene_rna_activity_log2fc_by_celltype.tsv")

NETWORK_EDGES <- file.path(NETWORK_DIR, "cortex_1_erf_target_edges_scored.tsv")
DRIVER_RANKING <- file.path(NETWORK_DIR, "cortex_1_erf_driver_ranking.tsv")
WEIGHTED_OUT_DEGREE <- file.path(NETWORK_DIR, "cortex_1_weighted_out_degree_source_data.tsv")

# Plotting parameters
TARGET_CELLTYPE <- "Cortex_1"

# Use neutral color-blind-friendly palettes only if needed.
# Users may adjust colors during final figure assembly.
CLADE_COLORS <- c(ECA = "#D9A441", ECB = "#7B4FA3", A = "#D9A441", B = "#7B4FA3")
CONDITION_COLORS <- c(Control = "#4D4D4D", Salinity = "#2B8CBE", Ctrl = "#4D4D4D", NaCl = "#2B8CBE")

dir_create(OUT_DIR, recurse = TRUE)
dir_create(PLOT_DIR, recurse = TRUE)
dir_create(SOURCE_DATA_DIR, recurse = TRUE)

message("PROJECT_DIR: ", PROJECT_DIR)
message("OUT_DIR: ", OUT_DIR)
