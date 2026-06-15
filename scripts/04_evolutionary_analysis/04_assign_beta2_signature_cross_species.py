#!/usr/bin/env python3

04_assign_beta2_signature_cross_species.py

Assign cross-species AP2/ERF beta2 residue signatures from an aligned FASTA.

This script requires manually verified alignment columns corresponding to the
two diagnostic beta2 residues. Columns are 1-based alignment positions.

Input FASTA header format:
    >Species|GeneID

or:
    >GeneID

If species is not encoded in the header, provide a metadata table separately
or edit the output after running.

Example:
    python scripts/04_evolutionary_analysis/04_assign_beta2_signature_cross_species.py \
      --aligned-fasta results/04_evolutionary_analysis/input/beta2/cross_species_erf_ap2_domains_aligned.fa \
      --columns 25,26

import argparse
import os
from pathlib import Path


def read_fasta(path):
    seqs = {}
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if name is not None:
            seqs[name] = "".join(chunks)
    return seqs


def classify(sig):
    sig = sig.upper()
    if sig == "AA":
        return "AA"
    if sig == "WV":
        return "WV"
    return "neither"


def parse_header(header):
    if "|" in header:
        species, gene = header.split("|", 1)
        return species.replace("_", " "), gene
    return "unknown", header


def main():
    project = os.environ.get("PROJECT_DIR", os.getcwd())
    out_dir = os.environ.get("OUT_DIR", os.path.join(project, "results", "04_evolutionary_analysis"))
    default_fasta = os.environ.get("CROSS_SPECIES_DBD_FASTA", os.path.join(out_dir, "input", "beta2", "cross_species_erf_ap2_domains_aligned.fa"))

    parser = argparse.ArgumentParser()
    parser.add_argument("--aligned-fasta", default=default_fasta)
    parser.add_argument("--columns", required=True, help="Comma-separated 1-based alignment columns, e.g. 25,26")
    parser.add_argument("--out-table", default=os.path.join(out_dir, "input", "beta2", "cross_species_beta2_signature.tsv"))
    args = parser.parse_args()

    cols = [int(x) for x in args.columns.split(",")]
    if len(cols) != 2:
        raise SystemExit("--columns must contain exactly two comma-separated positions.")

    seqs = read_fasta(args.aligned_fasta)
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_table, "w") as out:
        out.write("species\tgeneID\tbeta2_col1\tbeta2_col2\tbeta2_signature\tbeta2_signature_class\n")
        for header, seq in sorted(seqs.items()):
            species, gene_id = parse_header(header)
            residues = []
            for col in cols:
                residues.append(seq[col - 1].upper() if 1 <= col <= len(seq) else "NA")
            sig = "".join(residues)
            out.write(f"{species}\t{gene_id}\t{residues[0]}\t{residues[1]}\t{sig}\t{classify(sig)}\n")

    print(f"Saved cross-species beta2 signature table: {args.out_table}")
    print("Verify species names and beta2 columns before using the output for final plots.")


if __name__ == "__main__":
    main()
