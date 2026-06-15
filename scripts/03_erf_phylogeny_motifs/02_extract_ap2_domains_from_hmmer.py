#!/usr/bin/env python3


02_extract_ap2_domains_from_hmmer.py

Extract AP2/ERF DNA-binding domain sequences from a protein FASTA file using
HMMER --domtblout coordinates.

Usage from repository root:
    python scripts/03_erf_phylogeny_motifs/02_extract_ap2_domains_from_hmmer.py

Environment variables are read from 00_config.sh when using the run-all script.
You can also pass explicit arguments:
    python 02_extract_ap2_domains_from_hmmer.py \
        --protein-fasta reference/sorghum_dapseq_tf_proteins.fa \
        --domtblout results/03_erf_phylogeny_motifs/sorghum_erf_hmmer.domtblout \
        --erf-metadata results/03_erf_phylogeny_motifs/erf_metadata_prepared.tsv \
        --out-fasta results/03_erf_phylogeny_motifs/sorghum_erf_ap2_domains.fa

import argparse
import os
import re
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


def write_fasta(records, path, width=70):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            seq = seq.upper()
            for i in range(0, len(seq), width):
                out.write(seq[i:i+width] + "\n")


def read_erf_gene_ids(path):
    gene_ids = set()
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        if "geneID" not in header:
            raise ValueError(f"Expected geneID column in {path}")
        idx = header.index("geneID")
        for line in handle:
            if line.strip():
                gene_ids.add(line.rstrip("\n").split("\t")[idx])
    return gene_ids


def parse_domtblout(path, domain_regex):
    hits = []
    pattern = re.compile(domain_regex, re.IGNORECASE)
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 23:
                continue

            # HMMER domtblout fields:
            # target name, accession, tlen, query name, accession, qlen,
            # full E-value, full score, full bias, domain #, of, c-Evalue,
            # i-Evalue, score, bias, hmm from, hmm to, ali from, ali to,
            # env from, env to, acc, description
            target_name = parts[0]      # HMM model/domain name
            target_acc = parts[1]
            query_name = parts[3]       # protein sequence ID
            ali_from = int(parts[17])
            ali_to = int(parts[18])
            i_evalue = float(parts[12])
            score = float(parts[13])

            searchable = f"{target_name} {target_acc}"
            if pattern.search(searchable):
                hits.append({
                    "protein_id": query_name,
                    "domain_name": target_name,
                    "domain_acc": target_acc,
                    "ali_from": min(ali_from, ali_to),
                    "ali_to": max(ali_from, ali_to),
                    "i_evalue": i_evalue,
                    "score": score,
                })
    return hits


def main():
    env_project = os.environ.get("PROJECT_DIR", os.getcwd())
    env_out = os.environ.get("OUT_DIR", os.path.join(env_project, "results", "03_erf_phylogeny_motifs"))

    parser = argparse.ArgumentParser()
    parser.add_argument("--protein-fasta", default=os.environ.get("PROTEIN_FASTA", os.path.join(env_project, "reference", "sorghum_dapseq_tf_proteins.fa")))
    parser.add_argument("--domtblout", default=os.path.join(env_out, "sorghum_erf_hmmer.domtblout"))
    parser.add_argument("--erf-metadata", default=os.path.join(env_out, "erf_metadata_prepared.tsv"))
    parser.add_argument("--out-fasta", default=os.path.join(env_out, "sorghum_erf_ap2_domains.fa"))
    parser.add_argument("--out-table", default=os.path.join(env_out, "sorghum_erf_ap2_domain_coordinates.tsv"))
    parser.add_argument("--domain-regex", default=os.environ.get("HMM_DOMAIN_REGEX", "AP2|PF00847|AP2/ERF"))
    parser.add_argument("--concatenate-multiple", action="store_true", default=True, help="Concatenate multiple AP2 domains per protein if present.")
    args = parser.parse_args()

    Path(args.out_fasta).parent.mkdir(parents=True, exist_ok=True)

    seqs = read_fasta(args.protein_fasta)
    erf_ids = read_erf_gene_ids(args.erf_metadata)
    hits = parse_domtblout(args.domtblout, args.domain_regex)

    # Keep only ERF proteins listed in metadata when identifiers match.
    # If FASTA headers contain isoform suffixes, this still works when the geneID is contained in protein_id.
    filtered = []
    for h in hits:
        matched_gene = None
        for gid in erf_ids:
            if h["protein_id"] == gid or gid in h["protein_id"]:
                matched_gene = gid
                break
        if matched_gene is not None:
            h["geneID"] = matched_gene
            filtered.append(h)

    if not filtered:
        raise SystemExit("No AP2/ERF domain hits matched ERF metadata. Check FASTA headers, HMMER output and domain regex.")

    # For each gene, keep AP2 domain hits sorted by coordinate.
    by_gene = {}
    for h in filtered:
        by_gene.setdefault(h["geneID"], []).append(h)

    records = []
    coord_rows = []
    for gene_id, gene_hits in sorted(by_gene.items()):
        gene_hits = sorted(gene_hits, key=lambda x: (x["ali_from"], x["ali_to"], -x["score"]))

        # Remove nearly identical duplicate hits by same coordinate.
        unique_hits = []
        seen = set()
        for h in gene_hits:
            key = (h["ali_from"], h["ali_to"])
            if key not in seen:
                unique_hits.append(h)
                seen.add(key)

        protein_id = unique_hits[0]["protein_id"]
        if protein_id not in seqs:
            # try to find a sequence header containing the protein ID or gene ID
            candidates = [k for k in seqs if k == gene_id or gene_id in k or protein_id in k]
            if not candidates:
                print(f"WARNING: no FASTA sequence found for {gene_id} / {protein_id}; skipping")
                continue
            protein_id = candidates[0]

        protein_seq = seqs[protein_id]

        domain_seqs = []
        for i, h in enumerate(unique_hits, start=1):
            start = max(1, h["ali_from"])
            end = min(len(protein_seq), h["ali_to"])
            domain_seq = protein_seq[start-1:end].upper()
            domain_seqs.append(domain_seq)
            coord_rows.append([
                gene_id, protein_id, i, h["domain_name"], h["domain_acc"],
                start, end, h["i_evalue"], h["score"], domain_seq
            ])

        concatenated = "".join(domain_seqs)
        records.append((gene_id, concatenated))

    write_fasta(records, args.out_fasta)

    with open(args.out_table, "w") as out:
        out.write("\t".join([
            "geneID", "protein_id", "domain_index", "domain_name", "domain_acc",
            "protein_start", "protein_end", "i_evalue", "score", "domain_sequence"
        ]) + "\n")
        for row in coord_rows:
            out.write("\t".join(map(str, row)) + "\n")

    print(f"Saved AP2/ERF DBD FASTA: {args.out_fasta}")
    print(f"Saved domain coordinate table: {args.out_table}")
    print(f"Extracted domains for {len(records)} ERF genes")


if __name__ == "__main__":
    main()
