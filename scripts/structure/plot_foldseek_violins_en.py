"""
Figure for Section 4.5: structural similarity within vs between orthogroups
(15 AlphaFold3 models, 105 unordered pairs: 15 within, 90 between).

Reads the tables written by foldseek_statistics_unordered_en.py:
    results_tables/foldseek_unordered_pairs.csv
    results_tables/foldseek_unordered_summary.csv
Writes:
    figures/fig_structure_violins.png (300 dpi), .pdf and .svg

By default only panel A is drawn (one violin plot). With --full the
two supplementary panels B and C are added.

Panels:
    A  length-normalised TM-score (mean of query- and target-normalised),
       within vs between orthogroups; dashed line = 0.5 (same fold)
    B  sequence identity over the aligned region (Foldseek fident)
    C  TM-score against sequence identity, one point per pair

Usage (run from ~/project, environment with matplotlib, pandas, numpy):
    python3 scripts/structure/plot_foldseek_violins_en.py            # panel A only
    python3 scripts/structure/plot_foldseek_violins_en.py --full     # panels A, B, C
Colours: purple (within) and green (between), checked for colour-vision
deficiency (deuteranopia, protanopia, tritanopia) and contrast.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARGS = [a for a in sys.argv[1:] if not a.startswith("--")]
FULL = "--full" in sys.argv
PROJECT = Path(ARGS[0]) if ARGS else Path.home() / "project"
PAIRS = PROJECT / "results_tables" / "foldseek_unordered_pairs.csv"
SUMMARY = PROJECT / "results_tables" / "foldseek_unordered_summary.csv"
OUTDIR = PROJECT / "figures"
OUTDIR.mkdir(exist_ok=True)
NAME = "fig_structure_violins_full" if FULL else "fig_structure_violins"

C_WITHIN, C_BETWEEN = "#7B3294", "#1B9E77"   # purple / green, colour-blind-safe pair
rng = np.random.default_rng(1)

g = pd.read_csv(PAIRS)
S = pd.read_csv(SUMMARY).set_index("score")
w = g[g["type"] == "within"]
b = g[g["type"] == "between"]


def p_text(p):
    return "p < 10$^{-4}$" if p <= 1.01e-5 else f"p = {p:.2g}"


def violin(ax, col, ylabel, ylim, ref=None, title=None):
    data = [w[col].to_numpy(), b[col].to_numpy()]
    parts = ax.violinplot(data, positions=[1, 2], widths=0.8,
                          showextrema=False, showmedians=False)
    for body, c in zip(parts["bodies"], (C_WITHIN, C_BETWEEN)):
        body.set_facecolor(c); body.set_edgecolor(c)
        body.set_alpha(0.25); body.set_linewidth(1.2)
    for x, d, c in zip((1, 2), data, (C_WITHIN, C_BETWEEN)):
        ax.scatter(x + rng.uniform(-0.12, 0.12, len(d)), d, s=14, color=c,
                   alpha=0.85, edgecolor="white", linewidth=0.4, zorder=3)
        m = np.median(d)
        ax.hlines(m, x - 0.28, x + 0.28, color="black", linewidth=2, zorder=4)
    if ref is not None:
        ax.axhline(ref, color="grey", linestyle="--", linewidth=1, zorder=1)
    ax.set_xticks([1, 2])
    ax.set_xticklabels([f"Within\northogroup\n(n = {len(w)})",
                        f"Between\northogroups\n(n = {len(b)})"])
    ax.set_ylabel(ylabel)
    ax.set_ylim(*ylim)
    ax.set_title(title, loc="left", fontweight="bold", fontsize=11)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


if FULL:
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.4),
                             gridspec_kw={"width_ratios": [1, 1, 1.25]})
else:
    fig, ax_a = plt.subplots(1, 1, figsize=(4.6, 4.6))
    axes = [ax_a]

# A
ax = axes[0]
violin(ax, "tm_mean_norm", "TM-score (length-normalised)", (0.3, 1.12), ref=0.5,
       title="A" if FULL else None)
ax.text(0.02, 0.52, "TM-score 0.5", transform=ax.get_yaxis_transform(),
        fontsize=8, color="grey", va="bottom")
pA = S.loc["tm_mean_norm", "p_permutation_labels"]
ax.text(1.5, 1.05, f"permutation {p_text(pA)}", ha="center", fontsize=9)

if FULL:
    # B
    ax = axes[1]
    violin(ax, "fident", "Sequence identity over aligned region", (0.0, 1.12), title="B")
    pB = S.loc["fident", "p_permutation_labels"]
    ax.text(1.5, 1.05, f"permutation {p_text(pB)}", ha="center", fontsize=9)

    # C
    ax = axes[2]
    ax.scatter(b["fident"], b["tm_mean_norm"], s=18, color=C_BETWEEN, alpha=0.8,
               edgecolor="white", linewidth=0.4, label="Between orthogroups", zorder=2)
    ax.scatter(w["fident"], w["tm_mean_norm"], s=26, color=C_WITHIN, alpha=0.95,
               edgecolor="white", linewidth=0.4, label="Within orthogroup", zorder=3)
    ax.axhline(0.5, color="grey", linestyle="--", linewidth=1, zorder=1)
    ax.set_xlabel("Sequence identity over aligned region")
    ax.set_ylabel("TM-score (length-normalised)")
    ax.set_xlim(0, 1.0); ax.set_ylim(0.3, 1.05)
    ax.set_title("C", loc="left", fontweight="bold", fontsize=11)
    ax.legend(frameon=False, fontsize=8, loc="lower right")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)


fig.tight_layout()
for ext in ("png", "pdf", "svg"):
    fig.savefig(OUTDIR / f"{NAME}.{ext}", dpi=300,
                bbox_inches="tight")
print("Saved:", *(str(OUTDIR / f"{NAME}.{e}") for e in ("png", "pdf", "svg")), sep="\n  ")
print(f"matplotlib {matplotlib.__version__}, numpy {np.__version__}, pandas {pd.__version__}")
