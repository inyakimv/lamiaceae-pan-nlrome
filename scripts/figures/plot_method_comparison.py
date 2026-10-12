#!/usr/bin/env python3
"""
Figure 2 (Section 4.2): cross-validation between NLR-Annotator and ANNEVO+HMMER.

  a  Per-species number of complete NLR-Annotator candidates vs complete-like
     ANNEVO+HMMER candidates (dashed line: 1:1).
  b  Distribution of locus length (log scale) for the two methods, pooled across
     the 46 genomes.

Also prints every statistic quoted in the text: Pearson r, Spearman rho, Pearson r
without the largest-repertoire outlier, Mann-Whitney U, Cliff's delta, medians.

Inputs
  --counts   nlr_method_comparison_detailed.csv
             (columns: genome, nlr_annotator_complete, annevo_hmmer_complete_like)
  --nlr-len  gene_lengths_nlr_annotator.csv   (column: length_bp)
  --annevo-len gene_lengths_annevo.csv        (column: length_bp)
Outputs  <out>.png, <out>.pdf, <out>_source_data_a.csv, <out>_source_data_b_summary.csv

Usage
  python3 plot_method_comparison.py \
      --counts baseline_unfiltered/nlr_method_comparison_detailed.csv \
      --nlr-len gene_lengths_nlr_annotator.csv \
      --annevo-len gene_lengths_annevo.csv \
      --out Fig2_method_comparison
"""
import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

C_NLR = "#762A83"      # NLR-Annotator (purple; ColorBrewer PRGn, colour-blind safe with the green)
C_ANN = "#7FBF7B"      # ANNEVO+HMMER (green; CVD deltaE 35 vs purple, checked)
C_PT = "#4A6B8A"
OUTLIER_KEY = "Mentha_x_piperita"


def clean(ax):
    ax.grid(False)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.tick_params(labelsize=9)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--counts", required=True)
    ap.add_argument("--nlr-len", required=True)
    ap.add_argument("--annevo-len", required=True)
    ap.add_argument("--out", default="Fig2_method_comparison")
    a = ap.parse_args()

    d = pd.read_csv(a.counts)
    x = d["nlr_annotator_complete"].to_numpy()
    y = d["annevo_hmmer_complete_like"].to_numpy()
    r, pr = stats.pearsonr(x, y)
    rho, ps = stats.spearmanr(x, y)
    keep = ~d["genome"].str.contains(OUTLIER_KEY).to_numpy()
    r2, p2 = stats.pearsonr(x[keep], y[keep])

    ln = pd.read_csv(a.nlr_len)["length_bp"].to_numpy()
    la = pd.read_csv(a.annevo_len)["length_bp"].to_numpy()
    u, pu = stats.mannwhitneyu(ln, la, alternative="two-sided")
    delta = 2 * u / (len(ln) * len(la)) - 1  # Cliff's delta (NLR-Annotator vs ANNEVO)

    print("n species = %d | totals: NLR-Annotator %d, ANNEVO+HMMER %d" % (len(d), x.sum(), y.sum()))
    print("Pearson r = %.3f (p = %.2g); Spearman rho = %.3f (p = %.2g)" % (r, pr, rho, ps))
    print("Pearson r without %s = %.3f (p = %.2g)" % (OUTLIER_KEY, r2, p2))
    print("Loci: NLR-Annotator n = %d, median %.0f bp; ANNEVO+HMMER n = %d, median %.0f bp"
          % (len(ln), np.median(ln), len(la), np.median(la)))
    print("Mann-Whitney U = %.0f, p = %.3g; Cliff's delta = %.3f" % (u, pu, delta))

    fig, (axa, axb) = plt.subplots(1, 2, figsize=(11, 4.8), gridspec_kw={"width_ratios": [1, 1]})
    fig.patch.set_facecolor("white")

    # ---- panel a
    lim = max(x.max(), y.max()) * 1.05
    axa.plot([0, lim], [0, lim], ls="--", lw=1, color="#777777", zorder=1)
    axa.scatter(x, y, s=30, color=C_PT, edgecolor="white", linewidth=0.6, zorder=2)
    out = d[~keep].iloc[0]
    axa.annotate("$\\it{Mentha}$ ×\n$\\it{piperita}$",
                 (out["nlr_annotator_complete"], out["annevo_hmmer_complete_like"]),
                 xytext=(0, -9), textcoords="offset points", ha="center", va="top", fontsize=9, linespacing=1.15)
    axa.set_xlim(0, lim)
    axa.set_ylim(0, lim)
    axa.set_aspect("equal", adjustable="box")
    axa.set_xlabel("Complete NLR candidates, NLR-Annotator (loci)", fontsize=10)
    axa.set_ylabel("Complete-like NLR candidates, ANNEVO+HMMER (loci)", fontsize=10)
    axa.text(0.04, 0.96, "Pearson r = %.3f\nSpearman ρ = %.3f\nn = %d species" % (r, rho, len(d)),
             transform=axa.transAxes, va="top", ha="left", fontsize=9)
    clean(axa)

    # ---- panel b
    parts = axb.violinplot([ln, la], positions=[1, 2], showmedians=False, showextrema=False, widths=0.8)
    for body, col in zip(parts["bodies"], (C_NLR, C_ANN)):
        body.set_facecolor(col)
        body.set_edgecolor("#333333")
        body.set_linewidth(0.8)
        body.set_alpha(1.0)
    for pos, arr in ((1, ln), (2, la)):
        med = np.median(arr)
        axb.hlines(med, pos - 0.28, pos + 0.28, color="white", lw=2, zorder=3)
    axb.set_yscale("log")
    axb.set_xlim(0.4, 2.6)
    axb.set_xticks([1, 2])
    axb.set_xticklabels(["NLR-Annotator\n(n = %s)\nmedian %s bp" % (format(len(ln), ","), format(int(round(np.median(ln))), ",")),
                         "ANNEVO+HMMER\n(n = %s)\nmedian %s bp" % (format(len(la), ","), format(int(round(np.median(la))), ","))], fontsize=10)
    axb.set_ylabel("Locus length (bp, log scale)", fontsize=10)
    ymax = max(ln.max(), la.max())
    yb = ymax * 1.25
    axb.plot([1, 1, 2, 2], [yb / 1.12, yb, yb, yb / 1.12], color="#333333", lw=1)
    axb.text(1.5, yb * 1.06, "p < 0.001; Cliff's δ = %.2f" % delta, ha="center", va="bottom", fontsize=9)
    axb.set_ylim(top=yb * 2.0)
    clean(axb)

    for ax, lab in ((axa, "a"), (axb, "b")):
        ax.text(-0.12, 1.04, lab, transform=ax.transAxes, fontsize=14, fontweight="bold", va="bottom")

    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=300, facecolor="white")
    fig.savefig(a.out + ".pdf", facecolor="white")

    d[["genome", "nlr_annotator_complete", "annevo_hmmer_complete_like"]].to_csv(
        a.out + "_source_data_a.csv", index=False)
    pd.DataFrame({
        "method": ["NLR-Annotator", "ANNEVO+HMMER"],
        "n": [len(ln), len(la)],
        "median_bp": [np.median(ln), np.median(la)],
        "q1_bp": [np.percentile(ln, 25), np.percentile(la, 25)],
        "q3_bp": [np.percentile(ln, 75), np.percentile(la, 75)],
        "min_bp": [ln.min(), la.min()],
        "max_bp": [ln.max(), la.max()],
    }).to_csv(a.out + "_source_data_b_summary.csv", index=False)
    print("Saved", a.out + ".png/.pdf and source data")


if __name__ == "__main__":
    main()
