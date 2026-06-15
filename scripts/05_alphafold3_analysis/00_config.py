#!/usr/bin/env python3

"""Central Python configuration for AlphaFold3 analysis scripts."""

from pathlib import Path

PROJECT_DIR = Path.cwd()

METADATA_DIR = PROJECT_DIR / "metadata"
OUT_DIR = PROJECT_DIR / "results" / "05_alphafold3_analysis"
INPUT_DIR = OUT_DIR / "input"
AF3_OUTPUT_DIR = INPUT_DIR / "af3_outputs"

AF3_DESIGNS = METADATA_DIR / "af3_designs.tsv"
MOTIF_POSITION_MAP = METADATA_DIR / "af3_motif_position_map.tsv"
BETA2_RESIDUE_MAP = METADATA_DIR / "af3_beta2_residue_map.tsv"

AF3_JSON_DIR = OUT_DIR / "af3_input_json"
ATOM_TABLE_DIR = OUT_DIR / "af3_atom_tables"
PLOTS_DIR = OUT_DIR / "plots"

MANIFEST = OUT_DIR / "af3_output_manifest.tsv"
CONFIDENCE_SUMMARY = OUT_DIR / "af3_confidence_summary.tsv"

# Distance threshold, in Angstrom, used for protein-DNA heavy-atom contact counting.
CONTACT_DISTANCE_A = 4.0

# DNA atom classes used for contact summaries.
# Base atoms are used because the manuscript focuses on contacts to DNA bases at motif positions.
DNA_BACKBONE_ATOMS = {
    "P", "OP1", "OP2", "O1P", "O2P", "O3'", "O5'", "C5'", "C4'", "O4'", "C3'", "C2'", "C1'"
}

# Protein backbone atoms are excluded from beta2 side-chain distance summaries.
PROTEIN_BACKBONE_ATOMS = {"N", "CA", "C", "O", "OXT"}

for path in [OUT_DIR, INPUT_DIR, AF3_OUTPUT_DIR, AF3_JSON_DIR, ATOM_TABLE_DIR, PLOTS_DIR]:
    path.mkdir(parents=True, exist_ok=True)
