#!/usr/bin/env python3
"""
Species tree + horizontal bars of complete NLR candidates per species,
bars coloured by ploidy, plus (for numeric ploidy levels) a per-level
summary figure and statistics.

Two colour modes, both checked for colour-vision deficiency (CVD):

  --palette sequential   (default when all labels look like 2x, 4x, 6x ...)
        Ploidy is an ordered variable, so levels use ONE blue hue from light
        (low ploidy) to dark (high ploidy). Supports 2-5 levels (the ramps
        were validated for lightness steps; group levels if you have more).
        The level is also written in the species label, e.g. "Mentha aquatica (8x)",
        so colour is never the only cue.

  --palette categorical  (e.g. Diploid / Polyploid)
        Two classes: purple / orange (checked with the dataviz skill's
        validate_palette.js: CVD difference 32.6 on its scale; the orange is below 3:1
        against white, so the groups are also marked with a hatch unless --plain is used).
        --plain  draws only the coloured bars: no hatch, no legend (add --legend to keep it), no max-value label, no gridlines.
        With exactly two classes it also runs a Mann-Whitney U test
        and writes <out>_stats.txt and <out>_groups_summary.csv.
        Classes are ordered alphabetically (diploid, polyploid) unless --order is given.

Inputs
    --tree     Newick species tree. Tip names are matched to the tables after
               cleaning (GCF_/GCA_ prefix and "_genomic" removed, "_" -> space,
               case-insensitive).
    --counts   CSV, one row per genome (default columns: genome,
               nlr_annotator_complete  -> nlr_method_comparison_detailed.csv).
    --ploidy   CSV with columns: species, ploidy   (e.g. 2x, 4x, 6x, 8x)

Outputs (prefix --out)
    <out>.png / .pdf              tree + bars figure
    <out>_source_data.csv         exactly what was plotted
    <out>_levels.png / .pdf       counts per ploidy level (box + points)     [numeric levels]
    <out>_levels_summary.csv      n, median, quartiles, range per level      [numeric levels]
    stats printed to the terminal and saved in <out>_stats.txt               [numeric levels]

Usage
    python3 plot_nlr_ploidy.py --tree species_tree.nwk \
        --counts nlr_method_comparison_detailed.csv \
        --ploidy ploidy.csv --out nlr_counts_by_ploidy
"""
import argparse
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Patch
import numpy as np
import pandas as pd
from Bio import Phylo

# Categorical palette, slots 1-3 (validated for CVD separation):
CATEGORICAL = ["#6a3d9a", "#e08214"]   # purple (1st class), orange (2nd class)
HATCHES = ["", "////"]
# Ordinal blue ramps, light -> dark (validated: monotone lightness, visible steps,
# light end >= 2:1 against a white surface). Index = number of levels.
RAMPS = {
    2: ["#86b6ef", "#0d366b"],
    3: ["#86b6ef", "#2a78d6", "#0d366b"],
    4: ["#86b6ef", "#3987e5", "#1c5cab", "#0d366b"],
    5: ["#86b6ef", "#5598e7", "#2a78d6", "#184f95", "#0d366b"],
}
TEXT = "#0b0b0b"
MUTED = "#52514e"
GRID = "#e6e5e1"


def clean_name(raw):
    s = str(raw).strip()
    s = re.sub(r"^GC[AF]_\d+\.\d+_", "", s)
    s = re.sub(r"\.(faa|fasta|fa|pep|fna)$", "", s)
    s = re.sub(r"(\.annevo_proteins|_genomic)$", "", s)
    s = s.replace("_", " ")
    return re.sub(r"\s+", " ", s).strip()


def key(raw):
    return clean_name(raw).lower()


def level_value(label):
    """'4x' -> 4 ; anything else -> None."""
    m = re.match(r"^\s*(\d+)\s*x\s*$", str(label), re.I)
    return int(m.group(1)) if m else None


def natural_key(label):
    v = level_value(label)
    return (0, v, "") if v is not None else (1, 0, str(label).lower())


def layout_tree(tree):
    y = {leaf: i for i, leaf in enumerate(tree.get_terminals())}
    height = {}

    def walk(clade):
        if clade.is_terminal():
            height[clade] = 0
            return
        for c in clade.clades:
            walk(c)
        height[clade] = max(height[c] for c in clade.clades) + 1
        ys = [y[c] for c in clade.clades]
        y[clade] = (min(ys) + max(ys)) / 2

    walk(tree.root)
    return y, height


def draw_tree(ax, tree):
    y, height = layout_tree(tree)
    x = {c: -height[c] for c in height}

    def walk(clade):
        if clade.is_terminal():
            return
        ys = [y[c] for c in clade.clades]
        ax.plot([x[clade], x[clade]], [min(ys), max(ys)], color=MUTED, lw=1.0,
                solid_capstyle="butt")
        for c in clade.clades:
            ax.plot([x[clade], x[c]], [y[c], y[c]], color=MUTED, lw=1.0,
                    solid_capstyle="butt")
            walk(c)

    walk(tree.root)
    ax.set_xlim(x[tree.root] - 0.3, 0)
    ax.axis("off")


def style_axes(ax):
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)


def stats_block(data, classes):
    """Spearman (level vs count) and Kruskal-Wallis across levels."""
    from scipy import stats
    x = data["level"].to_numpy(dtype=float)
    y = data["count"].to_numpy(dtype=float)
    rho, p_rho = stats.spearmanr(x, y)
    groups = [data.loc[data["ploidy"] == c, "count"].to_numpy() for c in classes]
    try:
        h, p_kw = stats.kruskal(*groups)
    except ValueError:
        h, p_kw = float("nan"), float("nan")
    lines = [
        f"Species analysed: {len(data)}; ploidy levels: {', '.join(classes)}",
        f"Spearman correlation, ploidy level vs number of complete NLR candidates: "
        f"rho = {rho:.3f}, p = {p_rho:.3g}",
        f"Kruskal-Wallis test across ploidy levels: H = {h:.3f}, p = {p_kw:.3g}",
        "Caveat: species are not independent (shared phylogeny) and counts depend on how each",
        "assembly represents polyploid genomes (collapsed vs. phased subgenomes); treat these",
        "tests as descriptive evidence, not as formally independent-sample tests.",
    ]
    return rho, p_rho, h, p_kw, lines


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tree", required=True)
    ap.add_argument("--counts", required=True)
    ap.add_argument("--ploidy", required=True)
    ap.add_argument("--out", default="nlr_counts_by_ploidy")
    ap.add_argument("--plain", action="store_true",
                    help="only colour the bars: no hatch, no legend, no maximum-value label")
    ap.add_argument("--palette", choices=["auto", "sequential", "categorical"], default="auto")
    ap.add_argument("--genome-col", default="genome")
    ap.add_argument("--count-col", default="nlr_annotator_complete")
    ap.add_argument("--species-col", default="species")
    ap.add_argument("--ploidy-col", default="ploidy")
    ap.add_argument("--order", default=None,
                    help="class order, comma separated (categorical mode only)")
    ap.add_argument("--legend", action="store_true",
                    help="with --plain: still draw the legend (class names only, no counts)")
    ap.add_argument("--ladderize", action="store_true",
                    help="ladderize the tree (changes only the order of the tips)")
    ap.add_argument("--xlabel", default="Number of complete NLR candidates (gene count)")
    args = ap.parse_args()

    tree = Phylo.read(args.tree, "newick")
    if args.ladderize:
        tree.ladderize()
    tips = tree.get_terminals()

    counts = pd.read_csv(args.counts)
    counts["_key"] = counts[args.genome_col].map(key)
    count_by_key = dict(zip(counts["_key"], counts[args.count_col]))
    ploidy = pd.read_csv(args.ploidy)
    ploidy["_key"] = ploidy[args.species_col].map(key)
    ploidy_by_key = dict(zip(ploidy["_key"], ploidy[args.ploidy_col].astype(str).str.strip()))

    rows, miss_c, miss_p = [], [], []
    for t in tips:
        k, name = key(t.name), clean_name(t.name)
        if k not in count_by_key:
            miss_c.append(name)
        if k not in ploidy_by_key or ploidy_by_key[k] in ("", "nan"):
            miss_p.append(name)
        rows.append((name, count_by_key.get(k), ploidy_by_key.get(k)))
    if miss_c or miss_p:
        if miss_c:
            print("No count found for:", ", ".join(miss_c), file=sys.stderr)
        if miss_p:
            print("No ploidy found for:", ", ".join(miss_p), file=sys.stderr)
        sys.exit("Fix the tables (names must match the tree tips) and run again.")

    data = pd.DataFrame(rows, columns=["species", "count", "ploidy"])
    labels = list(data["ploidy"].unique())
    numeric = all(level_value(l) is not None for l in labels)
    mode = args.palette if args.palette != "auto" else ("sequential" if numeric else "categorical")
    if mode == "sequential" and not numeric:
        sys.exit("Sequential mode needs labels like 2x, 4x, 6x. Found: " + ", ".join(sorted(labels)))

    if mode == "sequential":
        classes = sorted(labels, key=natural_key)
        if len(classes) < 2 or len(classes) > max(RAMPS):
            sys.exit(f"{len(classes)} ploidy levels found; the validated ramps cover 2-{max(RAMPS)}. "
                     "Group the higher levels (e.g. 2x, 4x, 6x, >=8x labelled as 8x).")
        colour = dict(zip(classes, RAMPS[len(classes)]))
        hatch = {c: "" for c in classes}
        data["level"] = data["ploidy"].map(level_value)
        data["display"] = [f"{s} ({p})" for s, p in zip(data["species"], data["ploidy"])]
    else:
        classes = ([c.strip() for c in args.order.split(",")] if args.order
                   else sorted(labels, key=natural_key))
        if set(classes) != set(labels):
            sys.exit(f"--order {classes} does not match the labels {sorted(labels)}")
        if len(classes) > len(CATEGORICAL):
            sys.exit(f"{len(classes)} classes; categorical mode supports {len(CATEGORICAL)}. "
                     "Use ploidy labels like 2x, 4x, 6x for the sequential palette, or group.")
        colour = dict(zip(classes, CATEGORICAL))
        hatch = dict(zip(classes, HATCHES))
        if args.plain:
            hatch = {c: "" for c in classes}
        data["display"] = data["species"]

    # ---------------- figure 1: tree + bars ----------------
    n = len(data)
    plt.rcParams["hatch.linewidth"] = 1.3
    plt.rcParams["font.family"] = "DejaVu Sans"
    fig = plt.figure(figsize=(10.5, 0.25 * n + 1.2))
    gs = GridSpec(1, 3, width_ratios=[1.0, 1.15, 2.9], wspace=0.0, figure=fig,
                  left=0.01, right=0.985, top=0.99, bottom=0.07)
    ax_tree = fig.add_subplot(gs[0])
    ax_lab = fig.add_subplot(gs[1], sharey=ax_tree)
    ax_bar = fig.add_subplot(gs[2], sharey=ax_tree)

    draw_tree(ax_tree, tree)
    ax_tree.set_ylim(n - 0.5, -0.5)
    ax_lab.axis("off")
    for i, name in enumerate(data["display"]):
        ax_lab.text(0.97, i, name, ha="right", va="center", fontsize=9, color=TEXT)
    ax_lab.set_xlim(0, 1)

    ax_bar.barh(list(range(n)), data["count"], height=0.7,
                color=[colour[p] for p in data["ploidy"]],
                hatch=[hatch[p] for p in data["ploidy"]],
                edgecolor="white", linewidth=0, zorder=3)
    ax_bar.set_ylim(n - 0.5, -0.5)
    ax_bar.set_yticks([])
    style_axes(ax_bar)
    if not args.plain:
        ax_bar.xaxis.grid(True, color=GRID, lw=0.8, zorder=0)
        ax_bar.set_axisbelow(True)
    ax_bar.set_xlabel(args.xlabel, color=TEXT, fontsize=10)

    imax = int(data["count"].idxmax())
    if not args.plain:
      ax_bar.text(data["count"].iloc[imax] + data["count"].max() * 0.012, imax,
                f"{int(data['count'].iloc[imax])}", va="center", fontsize=9, color=TEXT)

    if (not args.plain) or args.legend:
        handles = [Patch(facecolor=colour[c], hatch=hatch[c], edgecolor="white", linewidth=0,
                         label=(c.capitalize() if args.plain else
                                f"{c} (n = {int((data['ploidy'] == c).sum())})")) for c in classes]
        leg = ax_bar.legend(handles=handles, title="Ploidy", loc="upper right", frameon=False,
                            fontsize=9, title_fontsize=9, handlelength=1.6, handleheight=1.1)
        for t in leg.get_texts():
            t.set_color(TEXT)
        leg.get_title().set_color(MUTED)

    fig.savefig(args.out + ".png", dpi=300)
    fig.savefig(args.out + ".pdf")
    plt.close(fig)
    data[["species", "ploidy", "count"]].rename(
        columns={"count": "complete_NLR_candidates"}).to_csv(args.out + "_source_data.csv", index=False)
    print(f"Wrote {args.out}.png, .pdf and _source_data.csv ({n} species; {mode} palette)")

    # ---------------- two-group comparison (categorical mode, 2 classes) ----------------
    if mode == "categorical" and len(classes) == 2:
        from scipy import stats as sps
        a_cls, b_cls = classes
        a = data.loc[data["ploidy"] == a_cls, "count"].to_numpy(dtype=float)
        b = data.loc[data["ploidy"] == b_cls, "count"].to_numpy(dtype=float)
        u, p_u = sps.mannwhitneyu(a, b, alternative="two-sided")
        summary = (data.groupby("ploidy")["count"]
                   .agg(n="size", median="median", q1=lambda s: s.quantile(.25),
                        q3=lambda s: s.quantile(.75), min="min", max="max")
                   .reindex(classes).reset_index())
        summary.to_csv(args.out + "_groups_summary.csv", index=False)
        lines = [f"{c}: n = {len(v)}, median = {np.median(v):.1f}, "
                 f"IQR = {np.percentile(v, 25):.1f}-{np.percentile(v, 75):.1f}, "
                 f"range = {v.min():.0f}-{v.max():.0f}" for c, v in ((a_cls, a), (b_cls, b))]
        lines.append(f"Mann-Whitney U = {u:.1f}, two-sided p = {p_u:.3g} (scipy {__import__('scipy').__version__})")
        text = "\n".join(lines)
        open(args.out + "_stats.txt", "w").write(text + "\n")
        print("\n" + text)
        print(f"\nWrote {args.out}_groups_summary.csv and {args.out}_stats.txt")

    # ---------------- figure 2 + stats (numeric levels) ----------------
    if mode == "sequential":
        rho, p_rho, h, p_kw, lines = stats_block(data, classes)
        summary = (data.groupby("ploidy")["count"]
                   .agg(n="size", median="median", q1=lambda s: s.quantile(.25),
                        q3=lambda s: s.quantile(.75), min="min", max="max")
                   .reindex(classes).reset_index())
        summary.to_csv(args.out + "_levels_summary.csv", index=False)

        fig2, ax = plt.subplots(figsize=(1.2 * len(classes) + 3.2, 5.0))
        rng = np.random.default_rng(1)
        for i, c in enumerate(classes):
            vals = data.loc[data["ploidy"] == c, "count"].to_numpy()
            bp = ax.boxplot([vals], positions=[i], widths=0.55, patch_artist=True,
                            showfliers=False,
                            medianprops=dict(color=TEXT, lw=1.6),
                            whiskerprops=dict(color=MUTED, lw=1.0),
                            capprops=dict(color=MUTED, lw=1.0),
                            boxprops=dict(facecolor=colour[c], edgecolor=MUTED, lw=1.0, alpha=0.55))
            ax.scatter(i + rng.uniform(-0.17, 0.17, len(vals)), vals, s=26,
                       color=colour[c], edgecolor=TEXT, linewidth=0.6, zorder=4)
        ax.set_xticks(range(len(classes)))
        ax.set_xticklabels([f"{c}\n(n = {int((data['ploidy'] == c).sum())})" for c in classes])
        ax.set_xlabel("Ploidy level", color=TEXT, fontsize=10)
        ax.set_ylabel(args.xlabel, color=TEXT, fontsize=10)
        style_axes(ax)
        ax.yaxis.grid(True, color=GRID, lw=0.8)
        ax.set_axisbelow(True)
        ax.text(0.02, 0.98, f"Spearman ρ = {rho:.2f}, p = {p_rho:.2g}\n"
                            f"Kruskal–Wallis p = {p_kw:.2g}",
                transform=ax.transAxes, va="top", fontsize=9, color=MUTED)
        fig2.tight_layout()
        fig2.savefig(args.out + "_levels.png", dpi=300)
        fig2.savefig(args.out + "_levels.pdf")
        plt.close(fig2)

        text = "\n".join(lines)
        open(args.out + "_stats.txt", "w").write(text + "\n")
        print("\n" + text)
        print(f"\nWrote {args.out}_levels.png, .pdf, _levels_summary.csv and _stats.txt")


if __name__ == "__main__":
    main()
