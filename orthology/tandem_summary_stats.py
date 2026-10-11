#!/usr/bin/env python3
"""
Summary of the tandem analysis (analyze_tandem_46species.py) and its correlation with the
number of complete NLR candidates per species (NLR-Annotator).
Usage: python3 tandem_summary_stats.py NLR_SUMMARY.csv TANDEM_ANALYSIS.csv
"""
import sys

import pandas as pd
from scipy import stats

nlr = pd.read_csv(sys.argv[1])
tan = pd.read_csv(sys.argv[2])
m = nlr.merge(tan, left_on="genoma", right_on="species", how="inner")
print("species:", len(m))
p = m["pct_genes_in_tandem_clusters"]
print("pct genes in tandem (top orthogroup): median %.1f, mean %.1f, min %.1f, max %.1f" % (p.median(), p.mean(), p.min(), p.max()))
print("copies of the top orthogroup per species: median %.1f, min %d, max %d" % (m.copies_in_species.median(), m.copies_in_species.min(), m.copies_in_species.max()))
r, pr = stats.pearsonr(m["NLR_completos"], p)
rho, ps = stats.spearmanr(m["NLR_completos"], p)
print("Pearson r = %.3f, p = %.3g; Spearman rho = %.3f, p = %.3g" % (r, pr, rho, ps))
print("species with tandem pct < 80:", int((p < 80).sum()), "| = 100:", int((p == 100).sum()))
