#!/usr/bin/env python3
"""Fig. 3 (Results 4.3): number of codon sites under selection per orthogroup.

Reads selection_summary.csv (written by selection_summary.py) and draws grouped
bars of the sites detected by FEL (pervasive positive selection, p <= 0.1) and
MEME (episodic diversifying selection, p <= 0.1), ordered from highest to lowest
total signal. Writes <out>.png, <out>.pdf and <out>_source_data.csv.

Usage: python3 plot_selection_summary.py --summary selection_summary.csv --out Fig3_selection
"""
import argparse, csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

C_FEL, C_MEME = "#D55E00", "#56B4E9"   # colour-blind-safe pair (Okabe-Ito vermillion / sky blue)
INK = "#222222"

# x-axis species labels (mathtext italics); edit here if the wording changes
LABELS = {
    "OG0000000": r"$\it{P.\ ornatus}$",
    "OG0000029": r"$\it{M.}$ × $\it{piperita}$",
    "OG0000018": r"$\it{M.}$ × $\it{piperita}$",
    "OG0000008": r"$\it{G.\ tetrahit}$",
    "OG0000038": r"$\it{S.\ hispanica}$",
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary", default="selection_summary.csv")
    ap.add_argument("--out", default="Fig3_selection")
    a = ap.parse_args()
    rows = list(csv.DictReader(open(a.summary)))
    rows.sort(key=lambda r: -(int(r["FEL_positive"]) + int(r["MEME"])))
    x = np.arange(len(rows)); w = 0.36
    fel = [int(r["FEL_positive"]) for r in rows]
    meme = [int(r["MEME"]) for r in rows]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "text.color": INK,
                         "axes.edgecolor": INK, "axes.labelcolor": INK,
                         "xtick.color": INK, "ytick.color": INK})
    fig, ax = plt.subplots(figsize=(8.6, 3.6))
    fig.patch.set_facecolor("white"); ax.set_facecolor("white")
    b1 = ax.bar(x - w / 2, fel, w, color=C_FEL, edgecolor=INK, linewidth=0.8, label="FEL (pervasive)")
    b2 = ax.bar(x + w / 2, meme, w, color=C_MEME, edgecolor=INK, linewidth=0.8, label="MEME (episodic)")
    for bars, vals in ((b1, fel), (b2, meme)):
        for b, v in zip(bars, vals):
            ax.text(b.get_x() + b.get_width() / 2, v + 4, str(v), ha="center", va="bottom", fontsize=8)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{r["orthogroup"]}\n{LABELS.get(r["orthogroup"], r["focal_species"])}' for r in rows])
    ax.set_ylabel("Number of selected sites")
    ax.set_ylim(0, max(meme) * 1.12)
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(False)
    fig.tight_layout()
    fig.savefig(a.out + ".png", dpi=300, facecolor="white")
    fig.savefig(a.out + ".pdf", facecolor="white")

    with open(a.out + "_source_data.csv", "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["orthogroup", "copies_analysed", "FEL_positive", "FEL_negative", "MEME", "shared_FEL_MEME"])
        for r in rows:
            cw.writerow([r["orthogroup"], r["copies_analysed"], r["FEL_positive"], r["FEL_negative"],
                         r["MEME"], r["shared_FEL_MEME"]])
    print("written", a.out + ".png/.pdf/_source_data.csv")


if __name__ == "__main__":
    main()
