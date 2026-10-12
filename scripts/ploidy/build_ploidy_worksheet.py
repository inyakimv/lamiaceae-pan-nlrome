#!/usr/bin/env python3
"""
Worksheet for assigning ploidy per assembly: assembly accession, NCBI link, CCDB counts,
and empty columns "ploidy" and "source" to fill in by hand (genome paper / BioProject).

Usage
  python3 build_ploidy_worksheet.py --genomes-dir 00_genomes --ccdb ccdb_summary.csv \
      --out ploidy_worksheet.csv
"""
import argparse
import glob
import os
import re

import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--genomes-dir", default="00_genomes")
ap.add_argument("--ccdb", required=True)
ap.add_argument("--out", default="ploidy_worksheet.csv")
a = ap.parse_args()

rows = []
for f in sorted(glob.glob(os.path.join(a.genomes_dir, "*.fna.fai"))):
    m = re.match(r"(GC[AF]_\d+\.\d+)_(.+)_genomic\.fna\.fai", os.path.basename(f))
    if m:
        acc, sp = m.group(1), m.group(2).replace("_", " ")
        rows.append((sp, acc, f"https://www.ncbi.nlm.nih.gov/datasets/genome/{acc}/"))
w = pd.DataFrame(rows, columns=["species", "assembly", "ncbi_url"])
c = pd.read_csv(a.ccdb)[["species", "mode_2n", "observed_2n"]]
w["key"] = w.species.str.lower().str.replace(" x ", " ")
c["key"] = c.species.str.lower().str.replace(" x ", " ")
w = w.merge(c.drop(columns="species"), on="key", how="left").drop(columns="key")
w["ploidy"] = ""
w["source"] = ""
w.to_csv(a.out, index=False)
print(len(w), "assemblies ->", a.out)
