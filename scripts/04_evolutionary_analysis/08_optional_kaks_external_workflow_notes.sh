#!/usr/bin/env bash

# 08_optional_kaks_external_workflow_notes.sh
# This is a documentation-only helper script showing the external Ka/Ks workflow.
# It is not executed by 99_run_all_evolutionary_analysis.sh.
#
# The exact implementation depends on the selected ortholog FASTA files,
# codon sequence naming and software installation. In the manuscript workflow,
# final statistical summaries are generated from a precomputed AP2/ERF DBD-level
# Ka/Ks table by scripts 03 and 04 in this folder.
#
# Typical external workflow:
#
# 1. For each sorghum–ortholog pair, extract protein DBD sequences.
# 2. Align protein DBDs, for example with ClustalW or MAFFT.
# 3. Extract corresponding CDS regions.
# 4. Convert protein alignment + CDS to codon alignment using PAL2NAL.
# 5. Estimate Ka, Ks and Ka/Ks using PAML/codeml or KaKs_Calculator.
# 6. Combine results into:
#
#    results/04_evolutionary_analysis/input/kaks/ap2_dbd_kaks_results.tsv
#
# with columns:
#
#    sorghum_geneID
#    ortholog_geneID
#    comparison
#    ka
#    ks
#    kaks
#
# Example skeleton commands:
#
# pal2nal.pl pair_protein.aln pair_cds.fa -output paml > pair.codon.phy
# codeml codeml.ctl
#
# or:
#
# KaKs_Calculator -i pair.axt -o pair.kaks.tsv -m NG
#
# Review and adapt this workflow to your local input files and installed software.
