#!/usr/bin/env python3

04_assign_beta2_signature_from_alignment.py

Assign the diagnostic beta2 residue signature from an aligned AP2/ERF DBD FASTA.

Important:
The exact alignment columns corresponding to the clade-defining beta2 residues
must be verified for each alignment. This script provides a reproducible way to
extract two positions once the columns are specified.

By default, the script attempts to use two alignment columns supplied through
--columns. Columns are 1-based alignment positions.

Example:
    python scripts/03_erf_phylogeny_motifs/04_assign_beta2_signature_from_alignment.py \
      --aligned-fasta results/03_erf_phylogeny_motifs/sorghum_erf_ap2_domains_aligned.fa \
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


def main():
    project = os.environ.get("PROJECT_DIR", os.getcwd())
    out_dir = os.environ.get("OUT_DIR", os.path.join(project, "results", "03_erf_phylogeny_motifs"))

    parser = argparse.ArgumentParser()
    parser.add_argument("--aligned-fasta", default=os.path.join(out_dir, "sorghum_erf_ap2_domains_aligned.fa"))
    parser.add_argument("--columns", default=None, help="Comma-separated 1-based alignment columns for the two beta2 residues, e.g. 25,26")
    parser.add_argument("--out-table", default=os.path.join(out_dir, "erf_beta2_signature.tsv"))
    args = parser.parse_args()

    if args.columns is None:
        raise SystemExit(
            "Please specify beta2 alignment columns using --columns, for example --columns 25,26. "
            "The columns should be verified from the AP2/ERF DBD alignment."
        )

    cols = [int(x) for x in args.columns.split(",")]
    if len(cols) != 2:
        raise SystemExit("--columns must contain exactly two comma-separated positions.")

    seqs = read_fasta(args.aligned_fasta)
    Path(args.out_table).parent.mkdir(parents=True, exist_ok=True)

    with open(args.out_table, "w") as out:
        out.write("geneID\tbeta2_col1\tbeta2_col2\tbeta2_signature\tbeta2_signature_class\n")
        for gene_id, seq in sorted(seqs.items()):
            residues = []
            for col in cols:
                if col < 1 or col > len(seq):
                    residues.append("NA")
                else:
                    residues.append(seq[col - 1].upper())
            sig = "".join(residues)
            out.write(f"{gene_id}\t{residues[0]}\t{residues[1]}\t{sig}\t{classify(sig)}\n")

    print(f"Saved beta2 signature table: {args.out_table}")
    print("Reminder: verify the selected columns against the curated ERF clade annotation.")


if __name__ == "__main__":
    main()
