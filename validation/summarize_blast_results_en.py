#!/usr/bin/env python3
"""
Reads the BLAST tabular output(s) against RefPlantNLR and
produces a single consolidated CSV summarizing, for each
query sequence, its best hit: which known/validated NLR gene
it resembles most, and at what % identity.

Usage:
    cd ~/project
    python3 summarize_blast_results_en.py

Automatically finds all *_blast_results.tsv files inside
10_validation/.
"""
import os
import glob
import csv
import re

PROJECT_DIR = os.path.expanduser("~/project")
OUT_DIR = os.path.join(PROJECT_DIR, "10_validation")

# BLAST outfmt 6 columns used:
# qseqid sseqid pident length evalue bitscore stitle
COLUMNS = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore", "stitle"]


def extract_species(qseqid):
    """qseqid format: OGxxxxxxx|GCA_..._1_Species_genomic_GCA_...1|locus|method"""
    parts = qseqid.split("|")
    if len(parts) < 2:
        return "unknown"
    m = re.search(r'_\d+_([A-Z][a-z]+_[a-z]+)_genomic', parts[1])
    return m.group(1).replace("_", " ") if m else "unknown"


def main():
    tsv_files = sorted(glob.glob(os.path.join(OUT_DIR, "*_blast_results.tsv")))
    if not tsv_files:
        print(f"ERROR: no *_blast_results.tsv files found inside {OUT_DIR}")
        return

    rows = []
    for tsv_path in tsv_files:
        og_id = os.path.basename(tsv_path).replace("_blast_results.tsv", "")
        with open(tsv_path) as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) != len(COLUMNS):
                    continue
                data = dict(zip(COLUMNS, cols))
                rows.append({
                    "orthogroup": og_id,
                    "species": extract_species(data["qseqid"]),
                    "query": data["qseqid"],
                    "best_hit_gene": data["stitle"],
                    "pident": float(data["pident"]),
                    "evalue": data["evalue"],
                    "bitscore": float(data["bitscore"]),
                })

    rows.sort(key=lambda r: (r["orthogroup"], -r["pident"]))

    out_path = os.path.join(OUT_DIR, "blast_validation_summary.csv")
    with open(out_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "orthogroup", "species", "query", "best_hit_gene", "pident", "evalue", "bitscore"
        ])
        writer.writeheader()
        writer.writerows(rows)

    print(f"Summary saved to: {out_path}")
    print(f"Total hits: {len(rows)}")

    # Quick stats per orthogroup
    by_og = {}
    for r in rows:
        by_og.setdefault(r["orthogroup"], []).append(r["pident"])
    print("\nPer orthogroup:")
    for og, pidents in by_og.items():
        avg = sum(pidents) / len(pidents)
        print(f"  {og}: {len(pidents)} hits, avg identity = {avg:.1f}%, "
              f"max identity = {max(pidents):.1f}%")


if __name__ == "__main__":
    main()
