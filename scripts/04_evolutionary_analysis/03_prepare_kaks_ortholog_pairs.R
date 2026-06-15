#!/usr/bin/env Rscript

# 03_prepare_kaks_ortholog_pairs.R
# Annotate precomputed AP2/ERF DBD Ka/Ks ortholog pairs with ERF clade labels.

source("scripts/04_evolutionary_analysis/00_config.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(stringr)
})

if (!file.exists(KAKS_RESULTS)) {
  message("KAKS_RESULTS not found: ", KAKS_RESULTS)
  message("Expected columns: sorghum_geneID, ortholog_geneID, comparison, ka, ks, kaks")
  message("Skipping Ka/Ks pair annotation.")
  quit(save = "no", status = 0)
}

if (!file.exists(ERF_CLADE_METADATA)) {
  stop("ERF clade metadata not found: ", ERF_CLADE_METADATA)
}

kaks <- read_tsv(KAKS_RESULTS, show_col_types = FALSE)
clade <- read_tsv(ERF_CLADE_METADATA, show_col_types = FALSE)

# Support common alternate column names.
if ("geneID" %in% colnames(kaks) && !("sorghum_geneID" %in% colnames(kaks))) {
  kaks <- kaks %>% rename(sorghum_geneID = geneID)
}
if ("KaKs" %in% colnames(kaks) && !("kaks" %in% colnames(kaks))) {
  kaks <- kaks %>% rename(kaks = KaKs)
}
if ("Ka" %in% colnames(kaks) && !("ka" %in% colnames(kaks))) {
  kaks <- kaks %>% rename(ka = Ka)
}
if ("Ks" %in% colnames(kaks) && !("ks" %in% colnames(kaks))) {
  kaks <- kaks %>% rename(ks = Ks)
}

required <- c("sorghum_geneID", "ortholog_geneID", "comparison", "ka", "ks", "kaks")
missing <- setdiff(required, colnames(kaks))
if (length(missing) > 0) {
  stop("Missing required columns in KAKS_RESULTS: ", paste(missing, collapse = ", "))
}

annotated <- kaks %>%
  mutate(
    ka = suppressWarnings(as.numeric(ka)),
    ks = suppressWarnings(as.numeric(ks)),
    kaks = suppressWarnings(as.numeric(kaks))
  ) %>%
  filter(is.finite(kaks), !is.na(sorghum_geneID), !is.na(ortholog_geneID)) %>%
  left_join(clade %>% select(geneID, clade), by = c("sorghum_geneID" = "geneID")) %>%
  mutate(
    clade_full = recode(clade, A = "ECA", B = "ECB", .default = NA_character_),
    valid_kaks = is.finite(kaks) & kaks >= 0
  ) %>%
  filter(valid_kaks, !is.na(clade)) %>%
  arrange(comparison, clade, sorghum_geneID, ortholog_geneID)

out_file <- file.path(OUT_DIR, "kaks_ortholog_pairs_annotated.tsv")
write_tsv(annotated, out_file)

message("Saved annotated Ka/Ks ortholog pairs: ", out_file)
message("Number of usable ortholog pairs: ", nrow(annotated))
print(table(annotated$comparison, annotated$clade, useNA = "ifany"))
message("Done.")
