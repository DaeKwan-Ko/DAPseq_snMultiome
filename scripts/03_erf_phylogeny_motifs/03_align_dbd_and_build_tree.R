#!/usr/bin/env Rscript

# 03_align_dbd_and_build_tree.R
# Align sorghum ERF AP2/ERF DNA-binding domains and build a neighbor-joining tree.

source("scripts/03_erf_phylogeny_motifs/00_config.R")

suppressPackageStartupMessages({
  library(msa)
  library(ape)
  library(Biostrings)
  library(readr)
  library(dplyr)
})

if (!file.exists(ERF_DBD_FASTA)) {
  stop("ERF DBD FASTA not found. Run 02_extract_ap2_domains_from_hmmer.py first: ", ERF_DBD_FASTA)
}

message("Reading DBD sequences: ", ERF_DBD_FASTA)
seqs <- readAAStringSet(ERF_DBD_FASTA)

if (length(seqs) < 3) {
  stop("At least 3 sequences are needed to build a tree.")
}

# ClustalW is used here to match the manuscript workflow.
# The msa package requires ClustalW to be available in the environment.
message("Aligning sequences with ClustalW through msa...")
aln <- msa(seqs, method = "ClustalW", order = "input")

aligned <- as(aln, "AAStringSet")
writeXStringSet(aligned, ERF_DBD_ALIGNED_FASTA)

message("Saved aligned FASTA: ", ERF_DBD_ALIGNED_FASTA)

# Convert alignment to DNAbin-like matrix for distance calculation.
aln_seqinr <- msaConvert(aln, type = "seqinr::alignment")
dist_mat <- dist.alignment(aln_seqinr, "identity")
tree <- nj(dist_mat)

# Preserve gene IDs as tip labels.
tree$tip.label <- names(seqs)[match(tree$tip.label, names(seqs), nomatch = seq_along(tree$tip.label))]

write.tree(tree, file = ERF_TREE_NWK)
message("Saved neighbor-joining tree: ", ERF_TREE_NWK)

# Save pairwise distance matrix for QC.
dist_df <- as.data.frame(as.matrix(dist_mat))
dist_df <- tibble::rownames_to_column(dist_df, "geneID")
write_tsv(dist_df, file.path(OUT_DIR, "sorghum_erf_dbd_pairwise_distance.tsv"))

message("Done.")
