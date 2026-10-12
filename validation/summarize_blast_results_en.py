#!/usr/bin/env python3
"""
Reads the BLAST tabular output(s) against RefPlantNLR and writes:

  10_validation/blast_validation_summary.csv
      one row per protein: its best hit (alignment with the highest bit
      score among the reported targets; ties broken by reference name),
      reference group, % identity, alignment length, query length and
      alignment coverage (alignment length / query length)
  10_validation/blast_validation_by_orthogroup.csv
      per orthogroup and overall: number of proteins, mean/max identity,
      mean coverage, and counts and percentages per reference group
      (source data of the BLAST table in the thesis)

Reference groups: NRC0; Rpi-amr1 (all allelic variants pooled); Rpi-chc1;
other NRC (other NRC helpers); other (any other RefPlantNLR protein).

Usage:
    cd ~/project
    python3 scripts/validation/summarize_blast_results_en.py

Expects BLAST outfmt 6 with the columns
qseqid sseqid pident length qlen slen evalue bitscore stitle
(written by run_blast_validation_en.sh).
"""
import os
import glob
import csv
import re
from collections import Counter

PROJECT_DIR = os.path.expanduser("~/project")
OUT_DIR = os.path.join(PROJECT_DIR, "10_validation")
COLUMNS = ["qseqid", "sseqid", "pident", "length", "qlen", "slen",
           "evalue", "bitscore", "stitle"]
GROUPS = ["NRC0", "Rpi-amr1", "Rpi-chc1", "other NRC", "other"]


def reference_group(gene):
    if gene.startswith("Rpi-amr1"):
        return "Rpi-amr1"
    if gene.startswith("Rpi-chc"):
        return "Rpi-chc1"
    if gene.startswith("NRC0"):
        return "NRC0"
    if "NRC" in gene:
        return "other NRC"
    return "other"


def extract_species(qseqid):
    """qseqid: OGxxxxxxx|GCA_..._1_Species_name_genomic_GCA_...1|locus|method"""
    parts = qseqid.split("|")
    if len(parts) < 2:
        return "unknown"
    m = re.match(r"^GC[AF]_\d+_\d_(.+?)_genomic_", parts[1])
    return m.group(1).replace("_", " ") if m else "unknown"


def summarise(rows):
    n = len(rows)
    out = {
        "proteins": n,
        "mean_pident": round(sum(r["pident"] for r in rows) / n, 1),
        "max_pident": round(max(r["pident"] for r in rows), 1),
        "mean_coverage": round(sum(r["coverage"] for r in rows) / n, 2),
        "coverage_ge_50pct": sum(1 for r in rows if r["coverage"] >= 0.5),
    }
    cnt = Counter(r["reference_group"] for r in rows)
    for g in GROUPS:
        out[g] = cnt.get(g, 0)
        out[g + "_pct"] = round(100 * cnt.get(g, 0) / n, 1)
    return out


def main():
    tsv_files = sorted(glob.glob(os.path.join(OUT_DIR, "*_blast_results.tsv")))
    if not tsv_files:
        raise SystemExit(f"ERROR: no *_blast_results.tsv files found inside {OUT_DIR}")

    best = {}
    for tsv_path in tsv_files:
        og_id = os.path.basename(tsv_path).replace("_blast_results.tsv", "")
        with open(tsv_path) as f:
            for line in f:
                if not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) != len(COLUMNS):
                    raise SystemExit(
                        f"ERROR: {tsv_path} has {len(cols)} columns, expected "
                        f"{len(COLUMNS)}. Re-run run_blast_validation_en.sh.")
                d = dict(zip(COLUMNS, cols))
                row = {
                    "orthogroup": og_id,
                    "species": extract_species(d["qseqid"]),
                    "query": d["qseqid"],
                    "best_hit_gene": d["stitle"],
                    "reference_group": reference_group(d["stitle"]),
                    "pident": float(d["pident"]),
                    "alignment_length": int(d["length"]),
                    "query_length": int(d["qlen"]),
                    "coverage": int(d["length"]) / int(d["qlen"]),
                    "evalue": d["evalue"],
                    "bitscore": float(d["bitscore"]),
                    "sseqid": d["sseqid"],
                }
                cur = best.get(row["query"])
                if (cur is None or row["bitscore"] > cur["bitscore"]
                        or (row["bitscore"] == cur["bitscore"]
                            and row["sseqid"] < cur["sseqid"])):
                    best[row["query"]] = row

    rows = sorted(best.values(), key=lambda r: (r["orthogroup"], r["query"]))

    out_path = os.path.join(OUT_DIR, "blast_validation_summary.csv")
    fields = ["orthogroup", "species", "query", "best_hit_gene", "reference_group",
              "pident", "alignment_length", "query_length", "coverage", "evalue", "bitscore"]
    with open(out_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        for r in rows:
            r = dict(r, coverage=round(r["coverage"], 3))
            w.writerow(r)
    print(f"Per-protein summary: {out_path} ({len(rows)} proteins)")

    by_og = {}
    for r in rows:
        by_og.setdefault(r["orthogroup"], []).append(r)
    og_fields = (["orthogroup", "proteins", "mean_pident", "max_pident",
                  "mean_coverage", "coverage_ge_50pct"]
                 + GROUPS + [g + "_pct" for g in GROUPS])
    og_path = os.path.join(OUT_DIR, "blast_validation_by_orthogroup.csv")
    with open(og_path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=og_fields)
        w.writeheader()
        for og in sorted(by_og):
            w.writerow(dict(orthogroup=og, **summarise(by_og[og])))
        w.writerow(dict(orthogroup="All", **summarise(rows)))
    print(f"Per-orthogroup summary: {og_path}")

    print("\nBy orthogroup (proteins, mean identity, NRC0 / Rpi-amr1 / Rpi-chc1 / other NRC / other):")
    for og in sorted(by_og):
        s = summarise(by_og[og])
        print(f"  {og}: {s['proteins']}, {s['mean_pident']}%, "
              + " / ".join(str(s[g]) for g in GROUPS))
    s = summarise(rows)
    print(f"  All: {s['proteins']}, {s['mean_pident']}%, "
          + " / ".join(str(s[g]) for g in GROUPS))
    print(f"  Coverage >= 50%: {s['coverage_ge_50pct']} of {s['proteins']}")

    print("\nMost frequent best-hit reference proteins:")
    for gene, n in Counter(r["best_hit_gene"] for r in rows).most_common(12):
        print(f"  {n:4d}  {gene}")


if __name__ == "__main__":
    main()
