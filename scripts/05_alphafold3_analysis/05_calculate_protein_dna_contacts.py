#!/usr/bin/env python3

"""
05_calculate_protein_dna_contacts.py

Calculate protein-DNA base contacts and minimum beta2-to-DNA-base distances.

Inputs:
    results/05_alphafold3_analysis/af3_atom_table_manifest.tsv
    metadata/af3_motif_position_map.tsv
    metadata/af3_beta2_residue_map.tsv

Outputs:
    protein_dna_base_contacts.tsv
    beta2_to_dna_base_min_distances.tsv
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


def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def read_atoms(path, cfg):
    atoms = []
    with open(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            element = row.get("element", "").upper()
            if element == "H":
                continue
            row["residue_number"] = int(row["residue_number"])
            row["x"] = float(row["x"])
            row["y"] = float(row["y"])
            row["z"] = float(row["z"])
            row["coord"] = np.array([row["x"], row["y"], row["z"]])
            atoms.append(row)
    return atoms


def dist(a, b):
    return float(np.linalg.norm(a["coord"] - b["coord"]))


def main():
    cfg = load_config()

    parser = argparse.ArgumentParser()
    parser.add_argument("--atom-table-manifest", default=str(cfg.OUT_DIR / "af3_atom_table_manifest.tsv"))
    parser.add_argument("--motif-position-map", default=str(cfg.MOTIF_POSITION_MAP))
    parser.add_argument("--beta2-residue-map", default=str(cfg.BETA2_RESIDUE_MAP))
    parser.add_argument("--contact-out", default=str(cfg.OUT_DIR / "protein_dna_base_contacts.tsv"))
    parser.add_argument("--beta2-out", default=str(cfg.OUT_DIR / "beta2_to_dna_base_min_distances.tsv"))
    args = parser.parse_args()

    manifest_path = Path(args.atom_table_manifest)
    motif_map_path = Path(args.motif_position_map)
    beta2_map_path = Path(args.beta2_residue_map)

    if not manifest_path.exists():
        raise SystemExit(f"Atom table manifest not found: {manifest_path}")

    if not motif_map_path.exists():
        raise SystemExit(
            f"Motif position map not found: {motif_map_path}. "
            "Create metadata/af3_motif_position_map.tsv."
        )

    if not beta2_map_path.exists():
        raise SystemExit(
            f"Beta2 residue map not found: {beta2_map_path}. "
            "Create metadata/af3_beta2_residue_map.tsv."
        )

    manifest = read_tsv(manifest_path)
    motif_map = read_tsv(motif_map_path)
    beta2_map = read_tsv(beta2_map_path)

    motif_by_design = {}
    for row in motif_map:
        motif_by_design.setdefault(row["design_id"], []).append({
            "dna_chain": row["dna_chain"],
            "dna_residue_number": int(row["dna_residue_number"]),
            "motif_position": row["motif_position"]
        })

    beta2_by_design = {}
    for row in beta2_map:
        beta2_by_design.setdefault(row["design_id"], []).append({
            "protein_chain": row["protein_chain"],
            "protein_residue_number": int(row["protein_residue_number"]),
            "beta2_site": row["beta2_site"]
        })

    contact_rows = []
    beta2_rows = []

    for row in manifest:
        design_id = row["design_id"]
        model_id = row["model_id"]

        if design_id not in motif_by_design:
            print(f"WARNING: no motif position map for {design_id}; skipping contacts.")
            continue

        atoms = read_atoms(row["atom_table"], cfg)

        protein_atoms = [
            a for a in atoms
            if a["residue_type"] == "protein"
            and a["atom_name"] not in cfg.PROTEIN_BACKBONE_ATOMS
        ]

        # DNA base atoms only, restricted to mapped motif positions.
        motif_positions = motif_by_design[design_id]
        dna_base_atoms = []
        for mp in motif_positions:
            for a in atoms:
                if (
                    a["residue_type"] == "DNA"
                    and a["chain_id"] == mp["dna_chain"]
                    and a["residue_number"] == mp["dna_residue_number"]
                    and a["atom_name"] not in cfg.DNA_BACKBONE_ATOMS
                ):
                    b = dict(a)
                    b["motif_position"] = mp["motif_position"]
                    dna_base_atoms.append(b)

        # Contact counts by motif position.
        for mp in motif_positions:
            target_atoms = [a for a in dna_base_atoms if a["motif_position"] == mp["motif_position"]]
            contact_pairs = []
            for pa in protein_atoms:
                for da in target_atoms:
                    d = dist(pa, da)
                    if d <= cfg.CONTACT_DISTANCE_A:
                        contact_pairs.append((pa, da, d))

            contacting_residues = sorted({
                f"{p['chain_id']}:{p['residue_name']}{p['residue_number']}"
                for p, _, _ in contact_pairs
            })

            contact_rows.append({
                "design_id": design_id,
                "model_id": model_id,
                "dna_chain": mp["dna_chain"],
                "dna_residue_number": mp["dna_residue_number"],
                "motif_position": mp["motif_position"],
                "contact_distance_threshold_A": cfg.CONTACT_DISTANCE_A,
                "n_protein_dna_base_atom_contacts": len(contact_pairs),
                "n_contacting_protein_residues": len(contacting_residues),
                "contacting_protein_residues": ";".join(contacting_residues)
            })

        # Minimum beta2 residue side-chain distance to DNA base atoms.
        if design_id in beta2_by_design:
            for bp in beta2_by_design[design_id]:
                beta_atoms = [
                    a for a in protein_atoms
                    if a["chain_id"] == bp["protein_chain"]
                    and a["residue_number"] == bp["protein_residue_number"]
                ]

                if not beta_atoms:
                    continue

                for mp in motif_positions:
                    target_atoms = [a for a in dna_base_atoms if a["motif_position"] == mp["motif_position"]]
                    min_d = math.nan
                    if target_atoms:
                        min_d = min(dist(pa, da) for pa in beta_atoms for da in target_atoms)

                    beta2_rows.append({
                        "design_id": design_id,
                        "model_id": model_id,
                        "protein_chain": bp["protein_chain"],
                        "protein_residue_number": bp["protein_residue_number"],
                        "beta2_site": bp["beta2_site"],
                        "dna_chain": mp["dna_chain"],
                        "dna_residue_number": mp["dna_residue_number"],
                        "motif_position": mp["motif_position"],
                        "min_beta2_to_dna_base_distance_A": min_d
                    })

    with open(args.contact_out, "w") as out:
        fieldnames = [
            "design_id", "model_id", "dna_chain", "dna_residue_number", "motif_position",
            "contact_distance_threshold_A", "n_protein_dna_base_atom_contacts",
            "n_contacting_protein_residues", "contacting_protein_residues"
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in contact_rows:
            writer.writerow(row)

    with open(args.beta2_out, "w") as out:
        fieldnames = [
            "design_id", "model_id", "protein_chain", "protein_residue_number",
            "beta2_site", "dna_chain", "dna_residue_number", "motif_position",
            "min_beta2_to_dna_base_distance_A"
        ]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in beta2_rows:
            writer.writerow(row)

    print(f"Saved protein-DNA base contacts: {args.contact_out}")
    print(f"Saved beta2-to-DNA base distances: {args.beta2_out}")


if __name__ == "__main__":
    main()
