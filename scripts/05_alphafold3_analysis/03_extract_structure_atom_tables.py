#!/usr/bin/env python3

"""
03_extract_structure_atom_tables.py

Parse AlphaFold3 mmCIF/CIF model structures into atom-level TSV tables.

Requires Biopython:
    pip install biopython

Outputs one atom table per model in:
    results/05_alphafold3_analysis/af3_atom_tables/

The atom table is used by downstream distance/RMSD scripts.
"""

import csv
import importlib.util
import re
from pathlib import Path


def load_config():
    cfg_path = Path(__file__).resolve().parent / "00_config.py"
    spec = importlib.util.spec_from_file_location("af3_config", cfg_path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)  # type: ignore
    return cfg


def require_biopython():
    try:
        from Bio.PDB import MMCIFParser  # noqa
        return MMCIFParser
    except ImportError as exc:
        raise SystemExit(
            "Biopython is required for parsing mmCIF files. Install with: pip install biopython"
        ) from exc


def sanitize_id(x):
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", str(x))


def is_dna_residue(resname):
    return resname.upper().strip() in {"DA", "DC", "DG", "DT", "A", "C", "G", "T"}


def residue_type(resname):
    return "DNA" if is_dna_residue(resname) else "protein"


def main():
    cfg = load_config()
    MMCIFParser = require_biopython()
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    import warnings

    if not cfg.MANIFEST.exists():
        raise SystemExit(f"Manifest not found. Run 02_collect_af3_outputs.py first: {cfg.MANIFEST}")

    parser = MMCIFParser(QUIET=True)

    rows_out = []
    with cfg.MANIFEST.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        manifest = [r for r in reader if r["file_type"] == "structure"]

    if not manifest:
        raise SystemExit("No structure rows found in manifest.")

    for row in manifest:
        design_id = row["design_id"]
        structure_path = Path(row["path"])
        model_id = "_".join(filter(None, [
            sanitize_id(design_id),
            f"seed{row.get('seed','')}" if row.get("seed") else "",
            f"sample{row.get('sample','')}" if row.get("sample") else "",
            f"model{row.get('model','')}" if row.get("model") else "",
            sanitize_id(structure_path.stem)
        ]))

        out_file = cfg.ATOM_TABLE_DIR / f"{model_id}.atoms.tsv"

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", PDBConstructionWarning)
            structure = parser.get_structure(model_id, str(structure_path))

        with out_file.open("w") as out:
            fieldnames = [
                "design_id", "model_id", "structure_file", "model_index",
                "chain_id", "residue_number", "insertion_code",
                "residue_name", "residue_type", "atom_name", "element",
                "x", "y", "z"
            ]
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()

            for model in structure:
                for chain in model:
                    for residue in chain:
                        hetflag, resseq, icode = residue.id
                        resname = residue.get_resname().strip()
                        for atom in residue:
                            writer.writerow({
                                "design_id": design_id,
                                "model_id": model_id,
                                "structure_file": str(structure_path),
                                "model_index": model.id,
                                "chain_id": chain.id,
                                "residue_number": resseq,
                                "insertion_code": icode.strip(),
                                "residue_name": resname,
                                "residue_type": residue_type(resname),
                                "atom_name": atom.get_name().strip(),
                                "element": (atom.element or "").strip(),
                                "x": f"{atom.coord[0]:.6f}",
                                "y": f"{atom.coord[1]:.6f}",
                                "z": f"{atom.coord[2]:.6f}"
                            })

        rows_out.append({
            "design_id": design_id,
            "model_id": model_id,
            "atom_table": str(out_file),
            "structure_file": str(structure_path)
        })

    map_file = cfg.OUT_DIR / "af3_atom_table_manifest.tsv"
    with map_file.open("w") as out:
        fieldnames = ["design_id", "model_id", "atom_table", "structure_file"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows_out:
            writer.writerow(row)

    print(f"Saved atom tables for {len(rows_out)} structures.")
    print(f"Saved atom table manifest: {map_file}")


if __name__ == "__main__":
    main()
