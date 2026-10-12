"""
Re-analysis of the Foldseek all-vs-all comparison of the 15 AlphaFold3
structures (5 orthogroups x 3 structures), counting every pair ONCE.

Why: Foldseek reports each pair in both directions (A->B and B->A), so
counting every row doubles n (30/180 instead of 15/90) and treats
non-independent rows as independent. Here:
  1. Each unordered pair is used once (scores of both directions averaged).
  2. Three scores are analysed: alntmscore (TM-score over the aligned
     region), and the length-normalised qtmscore/ttmscore (their mean and
     minimum per pair; the 0.5 same-fold threshold of Xu & Zhang 2010
     refers to length-normalised scores).
  3. Mann-Whitney U (within vs between) on unordered pairs, and a
     permutation test: the orthogroup labels are shuffled among the 15
     structures (3 per group kept), which respects that pairs share
     structures. Statistic = median(within) - median(between).
  4. A leave-one-orthogroup-out check and a per-structure summary.

Usage (run from ~/project, environment with pandas, numpy, scipy):
    python3 scripts/structure/foldseek_statistics_unordered_en.py
Output: results_tables/foldseek_unordered_*.csv and a printed summary.
"""
import re
import itertools
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu

PROJECT = Path.home() / "project"
TSV = PROJECT / "08_structure_comparison" / "structure_comparison.tsv"
OUT = PROJECT / "results_tables"
COLS = ["query", "target", "fident", "alnlen", "evalue", "bits",
        "qtmscore", "ttmscore", "alntmscore"]
N_PERM = 100000
SEED = 1


def og(name):
    m = re.match(r"(og\d+)_", name, re.IGNORECASE)
    return m.group(1).lower()


df = pd.read_csv(TSV, sep="\t", names=COLS)
df = df[df["query"] != df["target"]].copy()
df["pair"] = df.apply(lambda r: tuple(sorted((r["query"], r["target"]))), axis=1)
df["tm_mean_norm"] = (df["qtmscore"] + df["ttmscore"]) / 2
df["tm_min_norm"] = df[["qtmscore", "ttmscore"]].min(axis=1)

# average the two directions of every pair
g = df.groupby("pair").agg(
    n_rows=("query", "size"),
    fident=("fident", "mean"),
    alntmscore=("alntmscore", "mean"),
    tm_mean_norm=("tm_mean_norm", "mean"),
    tm_min_norm=("tm_min_norm", "mean"),
).reset_index()
g["a"] = g["pair"].apply(lambda p: p[0])
g["b"] = g["pair"].apply(lambda p: p[1])
g["og_a"] = g["a"].apply(og)
g["og_b"] = g["b"].apply(og)
g["type"] = np.where(g["og_a"] == g["og_b"], "within", "between")

structures = sorted(set(g["a"]) | set(g["b"]))
labels = np.array([og(s) for s in structures])
idx = {s: i for i, s in enumerate(structures)}
print(f"Structures: {len(structures)}; unordered pairs: {len(g)} "
      f"(within {int((g['type']=='within').sum())}, "
      f"between {int((g['type']=='between').sum())}); "
      f"rows per pair: {sorted(int(v) for v in g['n_rows'].unique())}")

ia = g["a"].map(idx).to_numpy()
ib = g["b"].map(idx).to_numpy()
rng = np.random.default_rng(SEED)
summary = []

for score in ["alntmscore", "tm_mean_norm", "tm_min_norm", "fident"]:
    x = g[score].to_numpy()
    w = g.loc[g["type"] == "within", score]
    b = g.loc[g["type"] == "between", score]
    u, p_mw = mannwhitneyu(w, b, alternative="two-sided")
    # rank-biserial (positive = within > between)
    rbc = 2 * u / (len(w) * len(b)) - 1
    obs = np.median(w) - np.median(b)
    count = 0
    for _ in range(N_PERM):
        perm = rng.permutation(labels)
        same = perm[ia] == perm[ib]
        d = np.median(x[same]) - np.median(x[~same])
        if abs(d) >= abs(obs) - 1e-12:
            count += 1
    p_perm = (count + 1) / (N_PERM + 1)
    summary.append({
        "score": score,
        "n_within": len(w), "n_between": len(b),
        "median_within": round(float(np.median(w)), 3),
        "median_between": round(float(np.median(b)), 3),
        "mean_within": round(float(w.mean()), 3),
        "mean_between": round(float(b.mean()), 3),
        "min_within": round(float(w.min()), 3),
        "min_between": round(float(b.min()), 3),
        "frac_ge_0.5_within": round(float((w >= 0.5).mean()), 3),
        "frac_ge_0.5_between": round(float((b >= 0.5).mean()), 3),
        "U": round(float(u), 1),
        "p_MannWhitney_unordered": float(f"{p_mw:.3g}"),
        "rank_biserial_within_minus_between": round(float(rbc), 3),
        "p_permutation_labels": float(f"{p_perm:.3g}"),
    })

S = pd.DataFrame(summary)
print("\nSUMMARY (every pair counted once)")
print(S.T.to_string(header=False))

# leave-one-orthogroup-out (drop all pairs involving that orthogroup)
rows = []
for score in ["tm_mean_norm", "alntmscore"]:
    for o in sorted(set(labels)):
        h = g[(g["og_a"] != o) & (g["og_b"] != o)]
        w = h.loc[h["type"] == "within", score]
        b = h.loc[h["type"] == "between", score]
        u, p = mannwhitneyu(w, b, alternative="two-sided")
        rows.append({"score": score, "dropped_orthogroup": o,
                     "n_within": len(w), "n_between": len(b),
                     "median_within": round(float(w.median()), 3),
                     "median_between": round(float(b.median()), 3),
                     "p_MannWhitney": float(f"{p:.3g}")})
L = pd.DataFrame(rows)
print("\nLEAVE-ONE-ORTHOGROUP-OUT")
print(L.to_string(index=False))

# per-orthogroup within-group scores
W = g[g["type"] == "within"][["og_a", "a", "b", "fident", "alntmscore",
                               "tm_mean_norm", "tm_min_norm"]].sort_values(["og_a", "alntmscore"])
print("\nWITHIN-ORTHOGROUP PAIRS")
print(W.round(3).to_string(index=False))

# per-structure mean score against the other 14
rows = []
for s in structures:
    m = g[(g["a"] == s) | (g["b"] == s)]
    rows.append({"structure": s, "orthogroup": og(s),
                 "mean_alntmscore_vs_others": round(float(m["alntmscore"].mean()), 3),
                 "mean_tm_mean_norm_vs_others": round(float(m["tm_mean_norm"].mean()), 3),
                 "mean_fident_vs_others": round(float(m["fident"].mean()), 3)})
P = pd.DataFrame(rows).sort_values("mean_alntmscore_vs_others")
print("\nPER-STRUCTURE MEANS (lowest first)")
print(P.to_string(index=False))

OUT.mkdir(exist_ok=True)
g.drop(columns=["pair"]).to_csv(OUT / "foldseek_unordered_pairs.csv", index=False)
S.to_csv(OUT / "foldseek_unordered_summary.csv", index=False)
L.to_csv(OUT / "foldseek_unordered_leave_one_og_out.csv", index=False)
P.to_csv(OUT / "foldseek_unordered_per_structure.csv", index=False)
print("\nSaved: results_tables/foldseek_unordered_{pairs,summary,leave_one_og_out,per_structure}.csv")
