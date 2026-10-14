#!/usr/bin/env python3
"""Fig. 4 (Results 4.3): circular gene tree of one orthogroup, coloured by species.

Reads the IQ-TREE ML tree of an orthogroup (UFBoot values as node labels),
midpoint-roots it and draws a circular phylogram. Branches whose descendants all
belong to one species ("monophyletic blocks") take that species' colour; mixed
branches are grey. An outer ring shows the species of each tip. Six key species
are coloured and labelled; all other species are grey (a colour-blind-safe
palette has a limited number of distinguishable hues).

Also writes:
  <out>_species_blocks.csv  per species: copies, monophyletic blocks, block sizes,
                            UFBoot support of the largest block
  <out>_tips.csv            tip order, species and ring colour (to redraw elsewhere)

Usage: python3 plot_og_tree.py --tree 04_phylogeny/OG0000000/OG0000000.treefile --out Fig4_OG0000000_tree
Requires Biopython and matplotlib.
"""
import argparse, csv, math, re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from Bio import Phylo

NAME = re.compile(r"^(GC[AF]_\d+_\d)_(.+?)_genomic_")
GREY, INK = "#C4C4C4", "#222222"
# six colours that stay distinguishable for every pair under colour-vision deficiency
COLOURS = {
    "Salvia divinorum": "#7A3B8F",
    "Salvia hispanica": "#eda100",
    "Salvia splendens": "#56B4E9",
    "Salvia miltiorrhiza": "#D55E00",
    "Plectranthus ornatus": "#CC79A7",
    "Isodon lophanthoides": "#2a78d6",
}


def species_of(name):
    m = NAME.match(name or "")
    return m.group(2).replace("_", " ") if m else "unknown"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tree", required=True)
    ap.add_argument("--out", default="Fig4_tree")
    a = ap.parse_args()

    tree = Phylo.read(a.tree, "newick")
    tree.root_at_midpoint()
    tree.ladderize()
    tips = tree.get_terminals()
    n = len(tips)
    sp = {t: species_of(t.name) for t in tips}
    depth = tree.depths()
    if max(depth.values()) == 0:      # no branch lengths: use topology
        depth = tree.depths(unit_branch_lengths=True)
    gap = math.radians(18)
    ang = {t: (i + 0.5) * (2 * math.pi - gap) / n for i, t in enumerate(tips)}

    leaves_of, species_set = {}, {}
    def collect(c):
        if c.is_terminal():
            leaves_of[c] = [c]
        else:
            leaves_of[c] = [l for ch in c.clades for l in collect(ch)]
        species_set[c] = {sp[l] for l in leaves_of[c]}
        return leaves_of[c]
    collect(tree.root)

    def angle(c):
        if c in ang: return ang[c]
        ang[c] = float(np.mean([angle(ch) for ch in c.clades]))
        return ang[c]
    angle(tree.root)

    def colour_of(c):
        s = species_set[c]
        return COLOURS.get(next(iter(s)), GREY) if len(s) == 1 else GREY

    rmax = max(depth[t] for t in tips)
    fig = plt.figure(figsize=(8.2, 8.2), facecolor="white")
    ax = fig.add_subplot(111, projection="polar", facecolor="white")
    ax.set_theta_zero_location("N"); ax.set_theta_direction(-1)
    ax.set_axis_off()

    def draw(c):
        for ch in c.clades:
            col = colour_of(ch)
            ax.plot([ang[ch]] * 2, [depth[c], depth[ch]], color=col, lw=0.9, solid_capstyle="butt")
            draw(ch)
        if c.clades:
            a0 = min(ang[ch] for ch in c.clades); a1 = max(ang[ch] for ch in c.clades)
            th = np.linspace(a0, a1, max(2, int((a1 - a0) / 0.01)))
            ax.plot(th, [depth[c]] * len(th), color=colour_of(c), lw=0.9)
    draw(tree.root)

    # outer ring: species of each tip
    ring0, ring_h = rmax * 1.03, rmax * 0.05
    ax.bar([ang[t] for t in tips], [ring_h] * n, width=(2 * math.pi - gap) / n * 1.02, bottom=ring0,
           color=[COLOURS.get(sp[t], GREY) for t in tips], edgecolor="none")

    # monophyletic blocks (maximal clades with a single species)
    blocks = {}
    def find_blocks(c):
        if len(species_set[c]) == 1:
            blocks.setdefault(next(iter(species_set[c])), []).append(c); return
        for ch in c.clades: find_blocks(ch)
    find_blocks(tree.root)

    # direct labels for the coloured species, at the middle of their largest block
    for s, col in COLOURS.items():
        if s not in blocks: continue
        big = max(blocks[s], key=lambda c: len(leaves_of[c]))
        th = float(np.mean([ang[l] for l in leaves_of[big]]))
        deg = (math.degrees(th)) % 360
        right = deg < 180
        ax.text(th, ring0 + ring_h + rmax * 0.04, s, rotation=(90 - deg) if right else (270 - deg),
                rotation_mode="anchor", ha="left" if right else "right", va="center",
                fontsize=8.5, style="italic", color=INK)

    # scale bar in the gap
    L = next(v for v in (0.05, 0.1, 0.2, 0.25, 0.5, 1.0) if v >= rmax / 8)
    th0 = 2 * math.pi - gap / 2
    ax.plot([th0, th0], [ring0 - L, ring0], color=INK, lw=1.6)
    ax.text(th0, ring0 - L * 1.15, f"{L:g}", ha="center", va="top", fontsize=7.5, color=INK)
    ax.set_ylim(0, rmax * 1.35)

    counts = {}
    for t in tips: counts[sp[t]] = counts.get(sp[t], 0) + 1
    handles = [Patch(facecolor=COLOURS[s], label=f"{s} ({counts.get(s, 0)})") for s in COLOURS if s in counts]
    others = [s for s in counts if s not in COLOURS]
    handles.append(Patch(facecolor=GREY, label=f"Other species ({len(others)}; {sum(counts[s] for s in others)} copies)"))
    leg = fig.legend(handles=handles, loc="lower left", frameon=False, fontsize=8.5, bbox_to_anchor=(0.02, 0.02))
    for txt in leg.get_texts()[:-1]: txt.set_style("italic")
    fig.text(0.98, 0.03, f"{n} copies, {len(counts)} species\nscale: substitutions per site",
             ha="right", va="bottom", fontsize=8, color=INK)
    fig.savefig(a.out + ".png", dpi=300, facecolor="white", bbox_inches="tight")
    fig.savefig(a.out + ".pdf", facecolor="white", bbox_inches="tight")

    with open(a.out + "_species_blocks.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["species", "copies", "monophyletic_blocks", "block_sizes", "largest_block_UFBoot"])
        for s in sorted(counts, key=lambda x: -counts[x]):
            bl = sorted(blocks.get(s, []), key=lambda c: -len(leaves_of[c]))
            sup = bl[0].confidence if bl and not bl[0].is_terminal() else ""
            w.writerow([s, counts[s], len(bl), " ".join(str(len(leaves_of[c])) for c in bl), sup])
    with open(a.out + "_tips.csv", "w", newline="") as f:
        w = csv.writer(f); w.writerow(["order", "tip", "species", "ring_colour"])
        for i, t in enumerate(tips, 1):
            w.writerow([i, t.name, sp[t], COLOURS.get(sp[t], GREY)])
    print("written", a.out + ".png/.pdf/_species_blocks.csv/_tips.csv")


if __name__ == "__main__":
    main()
