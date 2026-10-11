#!/usr/bin/env python3
"""
Assembly size, N50 and number of chromosome-scale sequences per genome, from the
samtools/SeqKit .fna.fai index files (no FASTA is read).

Usage
  python3 assembly_chromosome_counts.py --genomes-dir 00_genomes \
      --ccdb ccdb_summary.csv --out assembly_chromosome_counts.csv
  --min-chr-mb   length threshold for a "chromosome-scale" sequence (default 10 Mb)
The haploid number n_haploid_ccdb = mode_2n / 2 is taken from the CCDB summary
(fetch_ploidy_ccdb.py) and is only a rough comparison, not a ploidy call.
"""
import argparse
import glob
import os
import re

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genomes-dir", default="00_genomes")
ap.add_argument("--ccdb", default=None, help="ccdb_summary.csv from fetch_ploidy_ccdb.py")
ap.add_argument("--min-chr-mb", type=float, default=10.0)
ap.add_argument("--out", default="assembly_chromosome_counts.csv")
a = ap.parse_args()

rows = []
for f in sorted(glob.glob(os.path.join(a.genomes_dir, "*.fna.fai"))):
    m = re.match(r"(GC[AF]_\d+\.\d+)_(.+)_genomic\.fna\.fai", os.path.basename(f))
    if not m:
        continue
    acc, sp = m.group(1), m.group(2).replace("_", " ")
    L = pd.read_csv(f, sep="\t", header=None, usecols=[0, 1])[1].sort_values(ascending=False).values
    tot = int(L.sum())
    cum = 0
    n50 = 0
    for l in L:
        cum += l
        if cum >= tot / 2:
            n50 = int(l)
            break
    rows.append(dict(species=sp, assembly=acc, n_seqs=len(L), total_Mb=round(tot / 1e6, 1),
                     longest_Mb=round(L[0] / 1e6, 1), N50_Mb=round(n50 / 1e6, 1),
                     seqs_ge_min=int((L >= a.min_chr_mb * 1e6).sum()),
                     seqs_ge_1Mb=int((L >= 1e6).sum())))
d = pd.DataFrame(rows).rename(columns={"seqs_ge_min": f"seqs_ge_{a.min_chr_mb:g}Mb"})
if a.ccdb:
    c = pd.read_csv(a.ccdb)[["species", "mode_2n", "observed_2n"]]
    d["key"] = d.species.str.lower().str.replace(" x ", " ")
    c["key"] = c.species.str.lower().str.replace(" x ", " ")
    d = d.merge(c.drop(columns="species"), on="key", how="left").drop(columns="key")
    d["n_haploid_ccdb"] = d["mode_2n"] / 2
d.to_csv(a.out, index=False)
print(len(d), "assemblies ->", a.out)
