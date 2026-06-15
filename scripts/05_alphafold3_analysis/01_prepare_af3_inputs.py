#!/usr/bin/env python3

"""
01_prepare_af3_inputs.py

Generate AlphaFold3 JSON input files from metadata/af3_designs.tsv.

Expected input columns:
    design_id
    protein_variant
    protein_geneID
    clade
    protein_sequence
    dna_probe_id
    dna_sequence
    dna_motif_class
    notes
"""

import csv
import importlib.util
import json
import re
import sys
from pathlib import Path


def load_config():
    cfg_path = Path(__file__).resolve().parent / "00_config.py"
    spec = importlib.util.spec_from_file_location("af3_config", cfg_path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)  # type: ignore
    return cfg


def clean_sequence(seq, alphabet):
    seq = re.sub(r"\s+", "", seq).upper()
    bad = sorted(set(seq) - set(alphabet))
    if bad:
        raise ValueError(f"Unexpected sequence characters {bad} in sequence: {seq[:30]}...")
    return seq


def main():
    cfg = load_config()

    if not cfg.AF3_DESIGNS.exists():
        raise SystemExit(
            f"Missing {cfg.AF3_DESIGNS}. Create metadata/af3_designs.tsv first."
        )

    seeds = [1, 2, 3]

    with cfg.AF3_DESIGNS.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = [
            "design_id", "protein_variant", "protein_geneID", "clade",
            "protein_sequence", "dna_probe_id", "dna_sequence", "dna_motif_class"
        ]
        missing = [x for x in required if x not in reader.fieldnames]
        if missing:
            raise SystemExit(f"Missing required columns in {cfg.AF3_DESIGNS}: {', '.join(missing)}")

        n = 0
        summary_rows = []
        for row in reader:
            design_id = row["design_id"].strip()
            if not design_id:
                continue

            protein_sequence = clean_sequence(
                row["protein_sequence"],
                alphabet="ACDEFGHIKLMNPQRSTVWY"
            )
            dna_sequence = clean_sequence(
                row["dna_sequence"],
                alphabet="ACGT"
            )

            af3_job = {
                "name": design_id,
                "modelSeeds": seeds,
                "sequences": [
                    {
                        "protein": {
                            "id": "A",
                            "sequence": protein_sequence
                        }
                    },
                    {
                        "dna": {
                            "id": "B",
                            "sequence": dna_sequence
                        }
                    }
                ],
                "dialect": "alphafold3",
                "version": 1
            }

            out_json = cfg.AF3_JSON_DIR / f"{design_id}.json"
            with out_json.open("w") as out:
                json.dump(af3_job, out, indent=2)

            summary_rows.append({
                "design_id": design_id,
                "protein_variant": row["protein_variant"],
                "protein_geneID": row["protein_geneID"],
                "clade": row["clade"],
                "dna_probe_id": row["dna_probe_id"],
                "dna_motif_class": row["dna_motif_class"],
                "protein_length": len(protein_sequence),
                "dna_length": len(dna_sequence),
                "json_file": str(out_json)
            })
            n += 1

    summary_file = cfg.OUT_DIR / "af3_input_json_summary.tsv"
    with summary_file.open("w") as out:
        fieldnames = list(summary_rows[0].keys()) if summary_rows else [
            "design_id", "protein_variant", "protein_geneID", "clade",
            "dna_probe_id", "dna_motif_class", "protein_length", "dna_length", "json_file"
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in summary_rows:
            writer.writerow(row)

    print(f"Generated {n} AlphaFold3 JSON input files in {cfg.AF3_JSON_DIR}")
    print(f"Saved summary: {summary_file}")


if __name__ == "__main__":
    main()
