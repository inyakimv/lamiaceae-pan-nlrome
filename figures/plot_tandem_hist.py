#!/usr/bin/env python3
"""
Fig. 1b: distribution across species of the percentage of loci of the most expanded NLR
orthogroup that are organized in tandem clusters (output of analyze_tandem_46species.py).

Usage
  python3 plot_tandem_hist.py --tandem tandem_analysis_46species.csv --out Fig1b_tandem_hist
Writes <out>.png, <out>.pdf and <out>_source_data.csv
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--tandem", required=True)
ap.add_argument("--col", default="pct_genes_in_tandem_clusters")
ap.add_argument("--out", default="Fig1b_tandem_hist")
ap.add_argument("--bin", type=float, default=10.0, help="bin width in percentage points")
ap.add_argument("--color", default="#4A6B8A")
ap.add_argument("--xlabel", default="Tandem-clustered loci (%)")
a = ap.parse_args()

d = pd.read_csv(a.tandem)
v = d[a.col].astype(float)
edges = np.arange(0, 100 + a.bin, a.bin)
counts, _ = np.histogram(v, bins=edges)          # last bin includes 100
pd.DataFrame({"bin_start": edges[:-1], "bin_end": edges[1:], "n_species": counts}).to_csv(
    a.out + "_source_data.csv", index=False)

plt.rcParams["font.family"] = "DejaVu Sans"
fig, ax = plt.subplots(figsize=(5.2, 3.6))
ax.bar(edges[:-1], counts, width=a.bin, align="edge", color=a.color, edgecolor="white", linewidth=0.8)
ax.set_xlim(0, 100)
ax.set_xticks(np.arange(0, 101, 20))
ax.set_xlabel(a.xlabel, fontsize=10)
ax.set_ylabel("Number of species", fontsize=10)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.grid(False)
fig.patch.set_facecolor("white")
fig.tight_layout()
fig.savefig(a.out + ".png", dpi=300)
fig.savefig(a.out + ".pdf")
print(f"n = {len(v)} species; median {v.median():.1f}; counts per bin: {dict(zip([int(e) for e in edges[:-1]], [int(c) for c in counts]))}")
