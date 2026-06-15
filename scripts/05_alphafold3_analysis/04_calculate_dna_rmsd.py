#!/usr/bin/env python3

"""
04_calculate_dna_rmsd.py

Calculate DNA docking RMSD across independent AlphaFold3 predictions.

For each design_id, the first structure in the atom-table manifest is used as
the reference unless --reference-model-table is provided. DNA RMSD is computed
using matching DNA heavy atoms after Kabsch alignment.

Outputs:
    results/05_alphafold3_analysis/dna_rmsd_by_model.tsv
"""

import argparse
import csv
import importlib.util
import math
from pathlib import Path

import numpy as np


def load_config():
    cfg_path = Path(__file__).resolve().parent / "00_config.py"
    spec = importlib.util.spec_from_file_location("af3_config", cfg_path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)  # type: ignore
    return cfg


def read_atoms(path):
    atoms = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            if row["residue_type"] != "DNA":
                continue
            element = row.get("element", "").upper()
            if element == "H":
                continue
            key = (
                row["chain_id"],
                int(row["residue_number"]),
                row["residue_name"],
                row["atom_name"]
            )
            atoms.append((key, np.array([float(row["x"]), float(row["y"]), float(row["z"])])))
    return dict(atoms)


def kabsch_rmsd(P, Q):
    # P and Q are N x 3 arrays with corresponding atom coordinates.
    Pc = P - P.mean(axis=0)
    Qc = Q - Q.mean(axis=0)
    C = np.dot(Pc.T, Qc)
    V, S, Wt = np.linalg.svd(C)
    d = np.sign(np.linalg.det(np.dot(V, Wt)))
    D = np.diag([1.0, 1.0, d])
    U = np.dot(np.dot(V, D), Wt)
    P_rot = np.dot(Pc, U)
    diff = P_rot - Qc
    return math.sqrt((diff * diff).sum() / P.shape[0])


def main():
    cfg = load_config()
    parser = argparse.ArgumentParser()
    parser.add_argument("--atom-table-manifest", default=str(cfg.OUT_DIR / "af3_atom_table_manifest.tsv"))
    parser.add_argument("--out", default=str(cfg.OUT_DIR / "dna_rmsd_by_model.tsv"))
    args = parser.parse_args()

    manifest_path = Path(args.atom_table_manifest)
    if not manifest_path.exists():
        raise SystemExit(f"Atom table manifest not found. Run 03_extract_structure_atom_tables.py first: {manifest_path}")

    with manifest_path.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        manifest = list(reader)

    by_design = {}
    for row in manifest:
        by_design.setdefault(row["design_id"], []).append(row)

    out_rows = []
    for design_id, rows in sorted(by_design.items()):
        rows = sorted(rows, key=lambda r: r["model_id"])
        ref = rows[0]
        ref_atoms = read_atoms(ref["atom_table"])

        for row in rows:
            mob_atoms = read_atoms(row["atom_table"])
            common = sorted(set(ref_atoms) & set(mob_atoms))
            if len(common) < 3:
                rmsd = float("nan")
                n_atoms = len(common)
            else:
                Q = np.vstack([ref_atoms[k] for k in common])
                P = np.vstack([mob_atoms[k] for k in common])
                rmsd = kabsch_rmsd(P, Q)
                n_atoms = len(common)

            out_rows.append({
                "design_id": design_id,
                "model_id": row["model_id"],
                "reference_model_id": ref["model_id"],
                "n_matched_dna_heavy_atoms": n_atoms,
                "dna_rmsd_A": rmsd,
                "atom_table": row["atom_table"],
                "reference_atom_table": ref["atom_table"]
            })

    with open(args.out, "w") as out:
        fieldnames = [
            "design_id", "model_id", "reference_model_id",
            "n_matched_dna_heavy_atoms", "dna_rmsd_A",
            "atom_table", "reference_atom_table"
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in out_rows:
            writer.writerow(row)

    print(f"Saved DNA RMSD table: {args.out}")


if __name__ == "__main__":
    main()
