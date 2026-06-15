#!/usr/bin/env python3

"""
02_collect_af3_outputs.py

Collect AlphaFold3 output model files and confidence JSON files into a manifest.

The script recursively scans:
    results/05_alphafold3_analysis/input/af3_outputs/

It searches for:
    *.cif, *.mmcif
    *.json

Because AlphaFold3 Server and local AlphaFold3 output folders can have slightly
different naming conventions, users should manually inspect the manifest before
running downstream structure parsing.
"""

import csv
import importlib.util
import json
import re
from pathlib import Path


def load_config():
    cfg_path = Path(__file__).resolve().parent / "00_config.py"
    spec = importlib.util.spec_from_file_location("af3_config", cfg_path)
    cfg = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(cfg)  # type: ignore
    return cfg


def infer_design_id(path, known_designs):
    path_text = str(path)
    for d in sorted(known_designs, key=len, reverse=True):
        if d in path_text:
            return d

    # Fallback: use nearest parent folder name.
    for parent in [path.parent, path.parent.parent]:
        if parent.name:
            return parent.name
    return path.stem


def infer_seed_or_model(path):
    text = str(path)
    seed = None
    sample = None

    m = re.search(r"seed[-_]?(\d+)", text, flags=re.IGNORECASE)
    if m:
        seed = m.group(1)

    m = re.search(r"sample[-_]?(\d+)", text, flags=re.IGNORECASE)
    if m:
        sample = m.group(1)

    m = re.search(r"model[-_]?(\d+)", text, flags=re.IGNORECASE)
    model = m.group(1) if m else None

    return seed, sample, model


def read_known_designs(cfg):
    designs = []
    if cfg.AF3_DESIGNS.exists():
        with cfg.AF3_DESIGNS.open() as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            if "design_id" in reader.fieldnames:
                designs = [row["design_id"] for row in reader if row.get("design_id")]
    return designs


def parse_confidence_json(path):
    try:
        with path.open() as handle:
            data = json.load(handle)
    except Exception:
        return {}

    keys_of_interest = [
        "ranking_score", "iptm", "ptm", "fraction_disordered",
        "has_clash", "num_recycles"
    ]

    out = {}
    for key in keys_of_interest:
        if key in data and not isinstance(data[key], (list, dict)):
            out[key] = data[key]

    # Local AF3 summary_confidences may include chain_pair_iptm etc.
    for key, value in data.items():
        if key in keys_of_interest:
            continue
        if isinstance(value, (int, float, str, bool)):
            out[key] = value

    return out


def main():
    cfg = load_config()

    known_designs = read_known_designs(cfg)

    if not cfg.AF3_OUTPUT_DIR.exists():
        raise SystemExit(f"AF3 output directory not found: {cfg.AF3_OUTPUT_DIR}")

    structure_files = list(cfg.AF3_OUTPUT_DIR.rglob("*.cif")) + list(cfg.AF3_OUTPUT_DIR.rglob("*.mmcif"))
    json_files = list(cfg.AF3_OUTPUT_DIR.rglob("*.json"))

    rows = []
    for path in sorted(structure_files):
        design_id = infer_design_id(path, known_designs)
        seed, sample, model = infer_seed_or_model(path)
        rows.append({
            "design_id": design_id,
            "file_type": "structure",
            "seed": seed or "",
            "sample": sample or "",
            "model": model or "",
            "path": str(path),
            "filename": path.name
        })

    confidence_rows = []
    for path in sorted(json_files):
        design_id = infer_design_id(path, known_designs)
        seed, sample, model = infer_seed_or_model(path)
        rows.append({
            "design_id": design_id,
            "file_type": "json",
            "seed": seed or "",
            "sample": sample or "",
            "model": model or "",
            "path": str(path),
            "filename": path.name
        })

        conf = parse_confidence_json(path)
        if conf:
            conf_row = {
                "design_id": design_id,
                "seed": seed or "",
                "sample": sample or "",
                "model": model or "",
                "json_file": str(path)
            }
            conf_row.update(conf)
            confidence_rows.append(conf_row)

    with cfg.MANIFEST.open("w") as out:
        fieldnames = ["design_id", "file_type", "seed", "sample", "model", "path", "filename"]
        writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)

    if confidence_rows:
        # Allow flexible JSON keys.
        fieldnames = sorted(set().union(*(r.keys() for r in confidence_rows)))
        with cfg.CONFIDENCE_SUMMARY.open("w") as out:
            writer = csv.DictWriter(out, fieldnames=fieldnames, delimiter="\t")
            writer.writeheader()
            for row in confidence_rows:
                writer.writerow(row)

    print(f"Found {len(structure_files)} structure files and {len(json_files)} JSON files.")
    print(f"Saved manifest: {cfg.MANIFEST}")
    if confidence_rows:
        print(f"Saved confidence summary: {cfg.CONFIDENCE_SUMMARY}")


if __name__ == "__main__":
    main()
