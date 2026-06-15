# 07_peak_to_gene_networks

This folder contains scripts for constructing cell type-resolved ERF-to-target candidate regulatory networks by integrating DAP-seq binding, snMultiome peak-to-gene links, motif support and salinity-responsive RNA/activity behavior.

The workflow is designed for the Cortex_1-centered network analysis used in the manuscript, but the scripts can be applied to other annotated root cell types when the required inputs are available.

## Main goals

1. Run or document Signac `LinkPeaks` peak-to-gene linking by cell type and condition.
2. Scan linked ATAC peaks for ECA/ECB motif support.
3. Integrate DAP-seq ERF binding with linked peaks and candidate target genes.
4. Filter for salinity-responsive target-gene behavior using RNA and/or ATAC-derived gene activity.
5. Score ERF-to-target edges.
6. Build ECA-associated and ECB-associated directed candidate regulatory networks.
7. Rank ERF candidate drivers within matched-clade networks.
8. Export node and edge tables for Cytoscape, igraph or figure plotting.

## Expected upstream outputs

From `scripts/02_dapseq_annotation/`:

```text
results/02_dapseq_annotation/dapseq_candidate_peak_gene_links.tsv
results/02_dapseq_annotation/erf_candidate_targets_long.tsv
results/02_dapseq_annotation/erf_target_clade_counts.tsv
results/02_dapseq_annotation/erf_clade_biased_recurrent_targets.tsv
```

From `scripts/06_snmultiome_processing/`:

```text
results/06_snmultiome_processing/snmultiome_annotated.rds
results/06_snmultiome_processing/snmultiome_with_gene_activity.rds
results/06_snmultiome_processing/rna_log2fc_by_celltype.tsv
results/06_snmultiome_processing/activity_log2fc_by_celltype.tsv
results/06_snmultiome_processing/erf_rna_activity_log2fc_by_celltype.tsv
results/06_snmultiome_processing/target_gene_rna_activity_log2fc_by_celltype.tsv
```

Required metadata:

```text
metadata/erf_clade_annotation.tsv
metadata/salinity_linked_erfs.tsv
```

Optional metadata:

```text
metadata/representative_erfs.tsv
metadata/representative_targets.tsv
```

## Key parameters

The default network settings are:

```text
target cell type: Cortex_1
minimum within-clade ERFs for recurrent target definition: 6
maximum other-clade ERFs for recurrent target definition: 3
minimum absolute target RNA log2FC: 0
minimum absolute target activity log2FC: 0
```

The default target-behavior filter keeps genes with detectable RNA/activity summaries and non-zero salinity-associated change. Adjust thresholds in `00_config.R` as needed.

## Outputs

Default output directory:

```text
results/07_peak_to_gene_networks/
```

Main outputs:

- `linked_peaks_by_celltype.tsv`
- `linked_peak_motif_support.tsv`
- `erf_peak_gene_candidate_edges.tsv`
- `cortex1_erf_target_edges_scored.tsv`
- `cortex1_erf_target_nodes.tsv`
- `cortex1_eca_network_edges.tsv`
- `cortex1_ecb_network_edges.tsv`
- `cortex1_erf_driver_ranking.tsv`
- `representative_erf_target_network_edges.tsv`
- `network_output_qc.tsv`

## Suggested run order

From the repository root:

```bash
Rscript scripts/07_peak_to_gene_networks/01_linkpeaks_by_celltype.R
Rscript scripts/07_peak_to_gene_networks/02_scan_linked_peaks_for_motifs.R
Rscript scripts/07_peak_to_gene_networks/03_integrate_dapseq_linked_peaks.R
Rscript scripts/07_peak_to_gene_networks/04_filter_salinity_responsive_targets.R
Rscript scripts/07_peak_to_gene_networks/05_score_erf_target_edges.R
Rscript scripts/07_peak_to_gene_networks/06_build_clade_networks.R
Rscript scripts/07_peak_to_gene_networks/07_rank_erf_drivers.R
Rscript scripts/07_peak_to_gene_networks/08_export_network_tables.R
Rscript scripts/07_peak_to_gene_networks/09_qc_check_network_outputs.R
```

or run:

```bash
bash scripts/07_peak_to_gene_networks/99_run_all_peak_to_gene_networks.sh
```

## Notes

The resulting networks are candidate regulatory networks. Directed ERF-to-target edges indicate support from DAP-seq binding and snMultiome peak-to-gene linking; they do not by themselves prove in vivo TF occupancy or direct transcriptional regulation.

Large RDS objects and fragment files should not be committed to GitHub.

