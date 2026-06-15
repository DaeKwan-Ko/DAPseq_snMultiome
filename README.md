# DAPseq_snMultiome

This repository contains analysis scripts, metadata and documentation associated with the manuscript:

**A sorghum TF–DNA interactome links ancient ERF clade divergence to cell type-resolved salinity responses**

The repository provides workflows for DAP-seq processing, peak annotation, motif discovery, ERF clade classification, evolutionary analysis, AlphaFold3-derived structural analysis, single-nucleus multiome (snMultiome) processing, peak-to-gene linking, ERF regulatory network construction and key source-data plotting.

Raw sequencing data and large intermediate files are not stored in this repository. Instead, scripts and metadata are provided to reproduce the main analyses from deposited sequencing data, processed DAP-seq peak files and derived analysis tables.

Final manuscript figure assembly was performed manually from the analysis outputs. Scripts for key quantitative plots and source-data generation are included where applicable.

## Repository structure

```text
DAPseq_snMultiome/
├── docs/
├── environment/
├── metadata/
├── scripts/
│   ├── 01_dapseq_processing/
│   ├── 02_dapseq_annotation/
│   ├── 03_erf_phylogeny_motifs/
│   ├── 04_evolutionary_analysis/
│   ├── 05_alphafold3_analysis/
│   ├── 06_snmultiome_processing/
│   ├── 07_peak_to_gene_networks/
│   └── 08_key_plotting/
├── .gitignore
├── LICENSE
└── README.md
```

## Analysis modules

### 01. DAP-seq processing

`scripts/01_dapseq_processing/`

Documents post-JGI DAP-seq processing steps, including BAM correlation analysis, peak BED preparation, merged peak genome coverage, ChIPseeker-compatible peak preparation and motif discovery from top DAP-seq peaks.

### 02. DAP-seq annotation

`scripts/02_dapseq_annotation/`

Annotates DAP-seq peaks to the *Sorghum bicolor* v3.1.1 genome using ChIPseeker and defines conservative DAP-seq-associated candidate target genes. Candidate targets are assigned from promoter-proximal, genic or downstream peaks, while distal intergenic peaks are excluded from GO-based target assignment.

### 03. ERF phylogeny and motifs

`scripts/03_erf_phylogeny_motifs/`

Processes the 49 DAP-seq-profiled ERFs, including ERF metadata preparation, AP2/ERF DNA-binding domain extraction, DNA-binding domain-based phylogeny, MEME motif parsing and integration of ERF clade, motif and peak information.

### 04. Evolutionary analysis

`scripts/04_evolutionary_analysis/`

Summarizes ERF evolutionary analyses, including ortholog comparisons, AP2/ERF DNA-binding domain Ka/Ks values and cross-lineage β2 residue signatures across representative green plant species.

### 05. AlphaFold3 analysis

`scripts/05_alphafold3_analysis/`

Documents AlphaFold3-derived ERF–DNA structural analysis, including AlphaFold3 input preparation, model-output collection, DNA RMSD calculation, protein–DNA base-contact counting and β2-residue distance analysis.

### 06. snMultiome processing

`scripts/06_snmultiome_processing/`

Documents Cell Ranger ARC and Seurat/Signac analysis of sorghum root snMultiome data, including joint RNA/ATAC quality control, WNN integration, marker-based cell-type annotation, ATAC-derived gene activity calculation and RNA/activity summaries by cell type and condition.

### 07. Peak-to-gene networks

`scripts/07_peak_to_gene_networks/`

Integrates DAP-seq binding, snMultiome peak-to-gene links, motif support and salinity-responsive RNA/activity behavior to infer candidate ERF-to-target regulatory networks. The default network analysis focuses on Cortex_1.

### 08. Key plotting

`scripts/08_key_plotting/`

Generates reviewer-facing quantitative plots and source-data tables from processed outputs. These scripts reproduce the quantitative summaries used for figure generation, although final figure layout was assembled manually.

## Metadata

The `metadata/` folder contains lightweight tables used by the analysis scripts, including:

* `tf_list_142.tsv`: list of the 142 DAP-seq-profiled transcription factors and TF families.
* `erf_clade_annotation.tsv`: ECA/ECB clade annotation for the 49 profiled ERFs.
* `salinity_linked_erfs.tsv`: salinity-linked ERFs used for focused clade-level analyses.
* `sample_metadata_dapseq.tsv`: per-TF metadata for the DAP-seq resource.
* `sample_metadata_snmultiome.tsv`: sample-level metadata for control and salinity-treated snMultiome libraries.
* `species_list_beta2_signature.tsv`: species-level metadata for cross-lineage β2 signature analysis.

## Data availability

### Raw DAP-seq sequencing data

Raw DAP-seq sequencing data are available through the JGI Genome Portal.

* DOI: `10.46936/10.25585/60008787`

### Processed DAP-seq peak files and annotations

Processed DAP-seq peak files and peak annotation tables for the 142 sorghum transcription factors are available through Zenodo.

* DOI: `10.5281/zenodo.19615921`

### Raw snMultiome sequencing data

Raw snMultiome sequencing data are available through the NCBI Sequence Read Archive.

* BioProject: `PRJNA1454104`

Associated SRA records:

| Sample | Modality | SRA accession | BioSample      | Library ID           |
| ------ | -------- | ------------- | -------------- | -------------------- |
| NaCl1  | ATAC     | `SRR38138762` | `SAMN57292536` | `Sorghum_NaCl1_ATAC` |
| Ctl1   | ATAC     | `SRR38138763` | `SAMN57292535` | `Sorghum_Ctl1_ATAC`  |
| NaCl1  | GEX      | `SRR38138764` | `SAMN57292536` | `Sorghum_NaCl1_GEX`  |
| Ctl1   | GEX      | `SRR38138765` | `SAMN57292535` | `Sorghum_Ctl1_GEX`   |

## Reference genome and annotation

Sequence analyses used the *Sorghum bicolor* BTx623 reference genome and annotation.

* Reference genome: *Sorghum bicolor* v3.1.1
* Gene annotation: `Sbicolor_454_v3.1.1.gene.gff3`

## Software environment

Software versions and package information are summarized in:

```text
environment/software_versions.txt
```

Major software used in the analyses includes BBTools/BBDuk, Bowtie2, MACS3, MEME Suite, deepTools, bedtools, samtools, Cell Ranger ARC, Seurat, Signac, ChIPseeker, HMMER, ClustalW, OrthoFinder and AlphaFold3.

## Running the workflows

Each analysis folder contains its own `README.md`, configuration file and workflow-specific scripts. In general, scripts are intended to be run from the repository root.

Example:

```bash
bash scripts/01_dapseq_processing/99_run_all_dapseq_processing.sh
Rscript scripts/02_dapseq_annotation/01_make_txdb.R
```

Before running any workflow, edit the corresponding `00_config.R`, `00_config.sh` or `00_config.py` file to point to local reference files, downloaded sequencing data or processed input tables.

## Large files not included

This repository does not store:

* raw FASTQ files
* BAM/BAI files
* fragment files
* Cell Ranger ARC output directories
* large Seurat or Signac RDS objects
* large intermediate genome-wide files

These files should be obtained from the public data repositories listed above or generated locally using the provided scripts.

## License

This repository is released under the MIT License. See `LICENSE` for details.
