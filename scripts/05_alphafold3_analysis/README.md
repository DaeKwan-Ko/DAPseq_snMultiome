# 05_alphafold3_analysis

This folder contains scripts for AlphaFold3-derived ERF–DNA structural analyses.

The workflow documents how protein–DNA complex models were organized, how AlphaFold3 input JSON files can be generated, and how predicted protein–DNA complexes can be parsed to quantify DNA docking geometry, protein–DNA base contacts and minimum distances between diagnostic AP2/ERF β2 residues and DNA base atoms.

## Main goals

1. Prepare AlphaFold3 input JSON files for ERF protein–DNA combinations.
2. Collect AlphaFold3 output files from downloaded prediction folders.
3. Parse model confidence files and model structure files.
4. Calculate DNA docking RMSD across independent predictions.
5. Quantify protein contacts to DNA base atoms at motif positions.
6. Calculate minimum heavy-atom distances between β2 residues and DNA base atoms.
7. Summarize source data for Extended Data Fig. 4-style plots.

## Expected input metadata

Edit `00_config.py` and `00_config.R` before running.

The main input metadata file is:

```text
metadata/af3_designs.tsv
```

Expected columns:

```text
design_id
protein_variant
protein_geneID
clade
protein_sequence
dna_probe_id
dna_sequence
dna_motif_class
notes
```

Example `design_id` values:

```text
ECA_WT_ECA_motif_seedset
ECB_WT_ECB_motif_seedset
ECA_beta2_WV_ECB_motif_seedset
ECB_beta2_AA_ECA_motif_seedset
```

The `protein_sequence` column should contain the full protein sequence or the protein domain sequence submitted to AlphaFold3. Reciprocal β2 mutants should be represented as separate rows with the already-mutated sequence.

Optional helper metadata files:

```text
metadata/af3_motif_position_map.tsv
metadata/af3_beta2_residue_map.tsv
```

`af3_motif_position_map.tsv` expected columns:

```text
design_id
dna_chain
dna_residue_number
motif_position
```

`af3_beta2_residue_map.tsv` expected columns:

```text
design_id
protein_chain
protein_residue_number
beta2_site
```

## Expected AlphaFold3 outputs

Place downloaded AlphaFold3 model folders here:

```text
results/05_alphafold3_analysis/input/af3_outputs/
```

The scripts recursively search for:

- `.cif` or `.mmcif` model structure files
- `.json` confidence/summary files

The scripts are written to tolerate both AlphaFold3 Server-style downloaded folders and local AlphaFold3-style output folders, but exact filenames may vary. The output manifest should be checked manually before downstream analyses.

## Outputs

Default output directory:

```text
results/05_alphafold3_analysis/
```

Main outputs:

- `af3_input_json/`
- `af3_output_manifest.tsv`
- `af3_confidence_summary.tsv`
- `af3_atom_tables/`
- `dna_rmsd_by_model.tsv`
- `protein_dna_base_contacts.tsv`
- `beta2_to_dna_base_min_distances.tsv`
- `af3_structural_metrics_summary.tsv`
- `af3_extended_data_fig4_source_data.tsv`

## Suggested run order

From the repository root:

```bash
python scripts/05_alphafold3_analysis/01_prepare_af3_inputs.py
python scripts/05_alphafold3_analysis/02_collect_af3_outputs.py
python scripts/05_alphafold3_analysis/03_extract_structure_atom_tables.py
python scripts/05_alphafold3_analysis/04_calculate_dna_rmsd.py
python scripts/05_alphafold3_analysis/05_calculate_protein_dna_contacts.py
Rscript scripts/05_alphafold3_analysis/06_summarize_af3_metrics.R
Rscript scripts/05_alphafold3_analysis/07_qc_check_af3_outputs.R
```

or run:

```bash
bash scripts/05_alphafold3_analysis/99_run_all_alphafold3_analysis.sh
```

## Notes

AlphaFold3 predictions were used to compare relative docking behavior across matched protein–DNA combinations. These structural analyses are hypothesis-generating and should be interpreted together with DAP-seq motif preference and evolutionary analyses.

The scripts quantify distances using predicted structures. They do not claim experimental protein–DNA structures or direct biochemical validation.

