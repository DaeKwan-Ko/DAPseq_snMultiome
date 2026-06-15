# 04_evolutionary_analysis

This folder contains scripts for evolutionary analyses of AP2/ERF DNA-binding domain divergence, β2 residue signatures and Ka/Ks comparisons among ECA and ECB ortholog pairs.

The scripts are designed to document and reproduce the main comparative analyses used in the manuscript. They assume that ERF clade assignments and AP2/ERF DNA-binding domain analyses from `scripts/03_erf_phylogeny_motifs/` are already available.

## Main goals

1. Prepare species-level metadata for representative green plant lineages.
2. Parse OrthoFinder outputs to identify sorghum ERF ortholog pairs.
3. Prepare ortholog pair tables for sorghum–rice, sorghum–Arabidopsis and sorghum–soybean comparisons.
4. Summarize Ka/Ks values across the AP2/ERF DNA-binding domain.
5. Classify ERF β2 signatures across representative species as AA, WV or neither.
6. Build cross-lineage β2 signature count tables for plotting and interpretation.
7. Generate a combined AP2/ERF DBD phylogeny input table for broad plant-lineage analyses.

## Expected inputs

Edit `00_config.R` and `00_config.sh` before running.

Required metadata:

- `metadata/erf_clade_annotation.tsv`
  - columns: `geneID`, `clade`
- `metadata/species_list_beta2_signature.tsv`
  - columns: `lineage_order`, `species`, `taxonomy`, `major_group`, `total_erfs_analyzed`

Expected upstream outputs from `03_erf_phylogeny_motifs`:

- `results/03_erf_phylogeny_motifs/erf_metadata_prepared.tsv`
- `results/03_erf_phylogeny_motifs/sorghum_erf_ap2_domains_aligned.fa`
- `results/03_erf_phylogeny_motifs/erf_beta2_signature.tsv`

Expected OrthoFinder inputs:

- Orthogroup and orthologue tables from OrthoFinder.
- Recommended location:
  - `results/04_evolutionary_analysis/input/orthofinder/`

Expected Ka/Ks input:

A table containing AP2/ERF DBD-level Ka/Ks values for ortholog pairs. The expected minimum columns are:

```text
sorghum_geneID
ortholog_geneID
comparison
ka
ks
kaks
```

where `comparison` is one of:

```text
sorghum_rice
sorghum_arabidopsis
sorghum_soybean
```

Recommended filename:

```text
results/04_evolutionary_analysis/input/kaks/ap2_dbd_kaks_results.tsv
```

Expected cross-species β2 input:

A table containing ERF β2 signatures across species. The expected minimum columns are:

```text
species
geneID
beta2_signature
beta2_signature_class
```

where `beta2_signature_class` is one of `AA`, `WV` or `neither`.

If this table is not available, `04_assign_beta2_signature_cross_species.py` can be used after manually verifying the alignment columns corresponding to the β2 residues.

## Outputs

Default output directory:

```text
results/04_evolutionary_analysis/
```

Main outputs:

- `species_metadata_prepared.tsv`
- `ortholog_pairs_long.tsv`
- `kaks_ortholog_pairs_annotated.tsv`
- `kaks_summary_by_clade.tsv`
- `kaks_wilcoxon_tests.tsv`
- `cross_species_beta2_signature_counts.tsv`
- `cross_species_beta2_signature_long.tsv`
- `cross_lineage_beta2_summary_for_plotting.tsv`
- `combined_erf_phylogeny_metadata.tsv`

## Suggested run order

From the repository root:

```bash
Rscript scripts/04_evolutionary_analysis/01_prepare_species_metadata.R
Rscript scripts/04_evolutionary_analysis/02_parse_orthofinder_orthologues.R
Rscript scripts/04_evolutionary_analysis/03_prepare_kaks_ortholog_pairs.R
Rscript scripts/04_evolutionary_analysis/04_summarize_kaks_by_clade.R
Rscript scripts/04_evolutionary_analysis/05_build_cross_species_beta2_summary.R
Rscript scripts/04_evolutionary_analysis/06_prepare_combined_erf_phylogeny_metadata.R
Rscript scripts/04_evolutionary_analysis/07_qc_check_evolutionary_outputs.R
```

or run:

```bash
bash scripts/04_evolutionary_analysis/99_run_all_evolutionary_analysis.sh
```

## Notes

The Ka/Ks script summarizes precomputed Ka/Ks values. It does not attempt to rerun codon alignment, PAL2NAL or PAML/KaKs_Calculator by default, because those steps depend heavily on the exact ortholog FASTA files and external software environment. The repository therefore records the reproducible parsing, filtering and statistical summary workflow used after Ka/Ks values were generated.

The β2-signature extraction script requires manually verified alignment columns. This is intentional because AP2/ERF DBD alignment columns can shift depending on the sequence set and alignment method.
