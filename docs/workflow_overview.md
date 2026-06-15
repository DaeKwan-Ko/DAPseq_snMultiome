# Workflow overview

This repository contains analysis scripts and metadata associated with the manuscript:

**A sorghum TF–DNA interactome links ancient ERF clade divergence to cell type-resolved salinity responses**

The workflow is organized into eight analysis modules.

## 01. DAP-seq processing

Folder:

```text
scripts/01_dapseq_processing/
```

This module documents post-JGI DAP-seq processing steps, including BAM correlation analysis, peak BED preparation, merged peak genome coverage, preparation of peak files for ChIPseeker and motif discovery from top DAP-seq peaks.

## 02. DAP-seq annotation

Folder:

```text
scripts/02_dapseq_annotation/
```

This module annotates DAP-seq peaks to the *Sorghum bicolor* v3.1.1 genome using ChIPseeker and defines conservative DAP-seq-associated candidate target genes. Candidate targets are assigned from promoter-proximal, genic or downstream peaks, while distal intergenic peaks are excluded from GO-based target assignment.

## 03. ERF phylogeny and motifs

Folder:

```text
scripts/03_erf_phylogeny_motifs/
```

This module focuses on the 49 DAP-seq-profiled ERFs. It prepares ERF metadata, extracts AP2/ERF DNA-binding domains, builds DNA-binding domain-based phylogenies, parses MEME motif outputs and integrates ERF clade, motif and peak information.

## 04. Evolutionary analysis

Folder:

```text
scripts/04_evolutionary_analysis/
```

This module analyzes ERF clade divergence across species. It summarizes AP2/ERF DNA-binding domain Ka/Ks values, ortholog comparisons and cross-lineage β2 residue signatures across representative green plant lineages.

## 05. AlphaFold3 analysis

Folder:

```text
scripts/05_alphafold3_analysis/
```

This module documents AlphaFold3-derived ERF–DNA structural analysis. It prepares AlphaFold3 input JSON files, collects model outputs and calculates DNA docking RMSD, protein–DNA base-contact counts and minimum distances between β2 residues and DNA base atoms.

## 06. snMultiome processing

Folder:

```text
scripts/06_snmultiome_processing/
```

This module documents Cell Ranger ARC processing and Seurat/Signac analysis of sorghum root single-nucleus multiome data. It includes joint RNA/ATAC quality control, WNN integration, marker-based cell-type annotation, gene activity calculation and RNA/activity summaries by cell type and condition.

## 07. Peak-to-gene networks

Folder:

```text
scripts/07_peak_to_gene_networks/
```

This module integrates DAP-seq binding, snMultiome peak-to-gene links, motif support and salinity-responsive RNA/activity behavior to infer candidate ERF-to-target regulatory networks. The default analysis focuses on Cortex_1.

## 08. Key plotting

Folder:

```text
scripts/08_key_plotting/
```

This module generates key quantitative plots and source-data tables from the processed outputs. Final manuscript figure layout was assembled manually, but these scripts provide reproducible quantitative summaries for reviewer inspection.

## Large files

Raw FASTQ files, BAM files, fragment files, Cell Ranger ARC output directories, large RDS objects and other large intermediate files are not stored in this repository. Public accession numbers and data sources are listed in `docs/data_accession_links.md`.
