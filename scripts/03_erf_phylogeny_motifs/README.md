# 03_erf_phylogeny_motifs

This folder contains scripts for ERF-focused phylogenetic and motif analyses used to define and visualize ECA and ECB clade divergence.

The scripts start from metadata tables, ERF protein sequences, HMMER domain scans and MEME motif outputs generated from DAP-seq peaks. They are designed to document the analysis workflow rather than store large raw or intermediate files.

## Main goals

1. Prepare ERF metadata for the 49 sorghum DAP-seq-profiled ERFs.
2. Extract AP2/ERF DNA-binding domain sequences from protein FASTA files using HMMER output.
3. Align DNA-binding domains and infer a neighbor-joining phylogeny.
4. Parse MEME motif outputs from DAP-seq peaks.
5. Integrate ERF clade assignment, peak counts and motif summaries.
6. Generate motif logo commands for representative motifs.
7. Produce tables used for downstream ECA/ECB analyses.

## Expected inputs

Edit `00_config.sh` and `00_config.R` before running.

Required metadata:

- `metadata/tf_list_142.tsv`
  - columns: `geneID`, `tf_family`
- `metadata/erf_clade_annotation.tsv`
  - columns: `geneID`, `clade`
- `metadata/salinity_linked_erfs.tsv`
  - columns: `geneID`, `clade`

Required sequence/domain inputs:

- Protein FASTA for the 49 sorghum ERFs or all 142 DAP-seq TFs.
- HMMER `--domtblout` output from `hmmscan` or `hmmsearch`.
- A local HMM database such as Pfam-A, if HMMER scanning is being rerun.

Required motif inputs:

- MEME output directories or `meme.txt` / `.meme` files generated from DAP-seq peak sequences.
- These are usually generated from the top 100 peaks or summit-centered peak sequences.

Optional inputs:

- DAP-seq peak count summary from `scripts/01_dapseq_processing/`.
- Top motif annotation table if motifs were manually curated.

## Outputs

Default output directory:

```text
results/03_erf_phylogeny_motifs/
```

Main outputs:

- `erf_metadata_prepared.tsv`
- `sorghum_erf_ap2_domains.fa`
- `sorghum_erf_ap2_domains_aligned.fa`
- `sorghum_erf_dbd_tree.nwk`
- `erf_beta2_signature.tsv`
- `meme_top_motif_summary.tsv`
- `erf_phylogeny_motif_annotation_table.tsv`
- `ceqlogo_commands.sh`

## Suggested run order

From the repository root:

```bash
bash scripts/03_erf_phylogeny_motifs/01_run_hmmer_ap2_scan.sh
python scripts/03_erf_phylogeny_motifs/02_extract_ap2_domains_from_hmmer.py
Rscript scripts/03_erf_phylogeny_motifs/03_align_dbd_and_build_tree.R
python scripts/03_erf_phylogeny_motifs/04_assign_beta2_signature_from_alignment.py
Rscript scripts/03_erf_phylogeny_motifs/05_parse_meme_top_motifs.R
Rscript scripts/03_erf_phylogeny_motifs/06_integrate_phylogeny_motif_metadata.R
bash scripts/03_erf_phylogeny_motifs/07_make_ceqlogo_commands.sh
Rscript scripts/03_erf_phylogeny_motifs/08_qc_check_erf_outputs.R
```

or run:

```bash
bash scripts/03_erf_phylogeny_motifs/99_run_all_erf_phylogeny_motifs.sh
```

## Notes

The β2 signature script is intentionally configurable. Because AP2/ERF domain boundaries and alignment columns can shift slightly depending on protein set and alignment method, users should verify the extracted β2 residues against the curated ERF clade annotation before using the table for final biological interpretation.

Final figure layout and manual polishing of motif logos were performed outside this repository.

