# 08_key_plotting

This folder contains scripts for generating key quantitative plots and source-data tables from the processed DAP-seq, evolutionary, snMultiome and ERF network analysis outputs.

These scripts are not intended to reproduce the final Illustrator/PowerPoint layout exactly. Instead, they generate reviewer-facing plots and clean source-data tables for the main quantitative panels.

## Main goals

1. Generate key DAP-seq/ERF motif and clade summary plots.
2. Generate Ka/Ks and cross-lineage β2 signature summary plots.
3. Generate snMultiome RNA/activity and ERF-network summary plots.
4. Generate source-data tables for main figure-style panels.
5. Provide QC summaries confirming which upstream tables are available.

## Expected upstream outputs

From `scripts/02_dapseq_annotation/`:

```text
results/02_dapseq_annotation/peak_annotation_summary_by_tf.tsv
results/02_dapseq_annotation/erf_target_clade_counts.tsv
results/02_dapseq_annotation/erf_clade_biased_recurrent_targets.tsv
```

From `scripts/03_erf_phylogeny_motifs/`:

```text
results/03_erf_phylogeny_motifs/erf_phylogeny_motif_annotation_table.tsv
results/03_erf_phylogeny_motifs/sorghum_erf_dbd_tree.nwk
results/03_erf_phylogeny_motifs/meme_top_motif_summary.tsv
```

From `scripts/04_evolutionary_analysis/`:

```text
results/04_evolutionary_analysis/kaks_ortholog_pairs_annotated.tsv
results/04_evolutionary_analysis/kaks_summary_by_clade.tsv
results/04_evolutionary_analysis/kaks_wilcoxon_tests.tsv
results/04_evolutionary_analysis/cross_lineage_beta2_summary_for_plotting.tsv
```

From `scripts/06_snmultiome_processing/`:

```text
results/06_snmultiome_processing/celltype_condition_nucleus_counts.tsv
results/06_snmultiome_processing/erf_rna_activity_log2fc_by_celltype.tsv
results/06_snmultiome_processing/target_gene_rna_activity_log2fc_by_celltype.tsv
```

From `scripts/07_peak_to_gene_networks/`:

```text
results/07_peak_to_gene_networks/cortex_1_erf_target_edges_scored.tsv
results/07_peak_to_gene_networks/cortex_1_erf_driver_ranking.tsv
results/07_peak_to_gene_networks/cortex_1_weighted_out_degree_source_data.tsv
```

## Outputs

Default output directory:

```text
results/08_key_plotting/
```

Main outputs:

```text
plots/
source_data/
plotting_input_qc.tsv
```

Representative plot outputs:

```text
plots/fig2_erf_clade_target_counts.pdf
plots/fig3_kaks_pooled_boxplot.pdf
plots/fig3_cross_lineage_beta2_stacked_bar.pdf
plots/fig4_celltype_condition_nucleus_counts.pdf
plots/fig4_erf_rna_log2fc_by_celltype.pdf
plots/fig4_target_rna_activity_cortex1.pdf
plots/fig4_erf_driver_ranking_cortex1.pdf
```

## Suggested run order

From the repository root:

```bash
Rscript scripts/08_key_plotting/01_plot_erf_dapseq_summaries.R
Rscript scripts/08_key_plotting/02_plot_evolutionary_summaries.R
Rscript scripts/08_key_plotting/03_plot_snmultiome_summaries.R
Rscript scripts/08_key_plotting/04_plot_network_summaries.R
Rscript scripts/08_key_plotting/05_prepare_source_data_tables.R
Rscript scripts/08_key_plotting/06_qc_check_plotting_inputs.R
```

or run:

```bash
bash scripts/08_key_plotting/99_run_all_key_plotting.sh
```

## Notes

Final manuscript figures were assembled manually from analysis outputs. These scripts provide reproducible quantitative plots and source-data exports, but final figure typography, panel labels and layout may differ from the submitted figure files.

