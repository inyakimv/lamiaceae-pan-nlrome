#!/usr/bin/env python3
"""
Table 2 (assembly QC summary), genome size vs NLR correlations and NLR-count outliers.
Inputs: --qc 01_qc/qc_summary.csv, --counts nlr_method_comparison_detailed.csv
Outputs: <out>_table2.csv and <out>_stats.txt
"""
import argparse
import sys

import numpy as np
import pandas as pd
import scipy
from scipy import stats

ap = argparse.ArgumentParser()
ap.add_argument("--qc", required=True)
ap.add_argument("--counts", required=True)
ap.add_argument("--count-col", default="nlr_annotator_complete")
ap.add_argument("--out", default="qc_table2")
a = ap.parse_args()

q = pd.read_csv(a.qc)
q["assembly"] = q["assembly_file"].str.extract(r"^(GC[AF]_\d+\.\d+)")
q["species"] = (q["assembly_file"].str.replace(r"^GC[AF]_\d+\.\d+_", "", regex=True)
                .str.replace("_genomic", "").str.replace("_", " "))
q["Genome size (Mb)"] = q["sum_len_bp"] / 1e6
c = pd.read_csv(a.counts)
c["assembly"] = c["genome"].str.extract(r"^(GC[AF]_\d+\.\d+)")
m = q.merge(c[["assembly", a.count_col]], on="assembly", how="inner")
if len(m) != len(q):
    sys.exit(f"Only {len(m)} of {len(q)} assemblies matched between the two tables")
m["NLR density"] = m[a.count_col] / m["Genome size (Mb)"]

metrics = {"Genome size (Mb)": "Genome size (Mb)", "N50": "N50 (bp)", "GC_percent": "GC content (%)",
           "busco_complete_pct": "BUSCO complete (%)", "busco_single_pct": "BUSCO single-copy (%)",
           "busco_duplicated_pct": "BUSCO duplicated (%)", "busco_fragmented_pct": "BUSCO fragmented (%)",
           "busco_missing_pct": "BUSCO missing (%)", "NLR density": "NLR density (complete candidates / Mb)"}
rows = []
for k, label in metrics.items():
    s = m[k]
    rows.append([label, s.mean(), s.median(), s.std(ddof=1), s.min(), s.max()])
t2 = pd.DataFrame(rows, columns=["Metric", "Mean", "Median", "SD (n-1)", "Min", "Max"])
t2.to_csv(a.out + "_table2.csv", index=False)

lines = [f"n genomes: {len(m)}", f"counts column: {a.count_col}"]
r, p = stats.pearsonr(m["Genome size (Mb)"], m[a.count_col])
lines.append(f"Pearson genome size vs NLR count: r = {r:.3f}, p = {p:.2g}")
r2, p2 = stats.pearsonr(m["Genome size (Mb)"], m["NLR density"])
lines.append(f"Pearson genome size vs NLR density: r = {r2:.3f}, p = {p2:.2g}")
x = m[a.count_col]
q1, q3 = x.quantile(.25), x.quantile(.75)
fence = q3 + 1.5 * (q3 - q1)
lines.append(f"Tukey upper fence on NLR counts: {fence:.0f}")
out = m[x > fence][["species", a.count_col, "busco_complete_pct"]]
lines.append("Outliers:\n" + out.to_string(index=False))
low = m[m.busco_complete_pct < 95][["species", "busco_complete_pct", a.count_col]]
lines.append("BUSCO complete < 95%:\n" + low.to_string(index=False))
d = m["busco_duplicated_pct"]
lines.append(f"Duplicated BUSCO: >=90%: {(d >= 90).sum()} | 40-90%: {((d >= 40) & (d < 90)).sum()} | <40%: {(d < 40).sum()}")
lines.append(f"Python {sys.version.split()[0]}, pandas {pd.__version__}, NumPy {np.__version__}, SciPy {scipy.__version__}")
text = "\n".join(lines)
open(a.out + "_stats.txt", "w").write(text + "\n")
print(t2.round(3).to_string(index=False))
print("\n" + text)
