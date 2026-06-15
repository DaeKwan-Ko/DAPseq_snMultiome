# 01_dapseq_processing

Scripts for post-processing the JGI-provided sorghum DAP-seq dataset used in the manuscript.

The workflow starts from JGI-provided processed files, including aligned BAM files and peak files, rather than raw FASTQ files. Raw sequencing data are deposited through the JGI Genome Portal as described in the manuscript.

## Main inputs

- JGI-provided DAP-seq BAM files
- JGI-provided DAP-seq peak files
- Sorghum BTx623 reference genome FASTA
- Sorghum v3.1.1 GFF3 annotation
- Metadata table containing the 142 profiled TFs and TF-family assignments

## Scripts

1. `00_config.sh`  
   Central configuration file. Edit paths before running any script.

2. `01_generate_bam_filelists.sh`  
   Creates BAM file lists and sample-label files for all TFs and for ERF, MYB and WRKY subsets.

3. `02_multibam_summary_all.sb`  
   Runs deepTools `multiBamSummary` and `plotCorrelation` for all DAP-seq BAM files.

4. `03_multibam_summary_by_family.sb`  
   Runs family-specific deepTools correlation analysis for ERF, MYB and WRKY BAM files.

5. `04_prepare_peak_beds_for_chipseeker.sh`  
   Converts JGI/MACS-style peak files into BED files for ChIPseeker. Duplicated peak intervals are collapsed by keeping the highest score.

6. `05_peak_count_and_genome_coverage.sb`  
   Merges all peak BED files with bedtools and calculates genome coverage.

7. `06_chipseeker_annotation.R`  
   Builds a sorghum TxDb object from GFF3 and annotates DAP-seq peaks with ChIPseeker.

8. `07_extract_top_peaks_for_meme.sh`  
   Extracts top-ranked peak sequences for motif discovery with MEME.

9. `08_run_meme_motif_discovery.sb`  
   Runs MEME motif discovery on top-ranked DAP-seq peak sequences.

## Notes

Final figure assembly and manual formatting are not included here. These scripts reproduce the main DAP-seq processing, peak annotation and motif-discovery steps used for downstream analysis.

