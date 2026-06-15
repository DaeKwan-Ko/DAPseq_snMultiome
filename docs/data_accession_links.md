# Data accession links

This repository accompanies the manuscript:

**A sorghum TF–DNA interactome links ancient ERF clade divergence to cell type-resolved salinity responses**

## Raw DAP-seq sequencing data

Raw DAP-seq sequencing data are available through the JGI Genome Portal.

* Repository: JGI Genome Portal
* DOI: `10.46936/10.25585/60008787`
* Dataset: DAP-seq profiling of abiotic stress-responsive *Sorghum bicolor* transcription factors

## Processed DAP-seq peak files and annotations

Processed DAP-seq peak files and peak annotation tables for the 142 sorghum transcription factors profiled in this study are available through Zenodo.

* Repository: Zenodo
* DOI: `10.5281/zenodo.19615921`

## Raw snMultiome sequencing data

Raw single-nucleus multiome sequencing data are available through the NCBI Sequence Read Archive.

* Repository: NCBI SRA
* BioProject: `PRJNA1454104`

### snMultiome SRA records

| Sample | Modality | SRA accession | BioSample      | Library ID           |
| ------ | -------- | ------------- | -------------- | -------------------- |
| NaCl1  | ATAC     | `SRR38138762` | `SAMN57292536` | `Sorghum_NaCl1_ATAC` |
| Ctl1   | ATAC     | `SRR38138763` | `SAMN57292535` | `Sorghum_Ctl1_ATAC`  |
| NaCl1  | GEX      | `SRR38138764` | `SAMN57292536` | `Sorghum_NaCl1_GEX`  |
| Ctl1   | GEX      | `SRR38138765` | `SAMN57292535` | `Sorghum_Ctl1_GEX`   |

## Reference genome and annotation

Sequence analyses used the *Sorghum bicolor* BTx623 reference genome and gene annotation.

* Reference genome: *Sorghum bicolor* v3.1.1
* Gene annotation: `Sbicolor_454_v3.1.1.gene.gff3`

## External comparative datasets

Arabidopsis DAP-seq datasets used for comparative ERF motif analyses were obtained from the Plant Cistrome Database.

## Repository contents

This GitHub repository stores analysis scripts, metadata and documentation. It does not store raw FASTQ files, BAM files, Cell Ranger ARC output directories, fragment files, large Seurat objects or other large intermediate files.
