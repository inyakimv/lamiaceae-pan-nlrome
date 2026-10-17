#!/usr/bin/env python3
"""Clade-level ML tree (IQ-TREE) figure: ADR1/NRG1 clade of Lamiaceae with reference NLRs.

Usage: python3 07_fig_clade_tree.py TREEFILE OUT_PREFIX [collapsed|full]
TREEFILE: IQ-TREE .treefile with node labels 'SH-aLRT/UFBoot' (e.g. Clade18_tree.treefile).
Leaves: <Species>|<gene> (Lamiaceae), REF_<name> (references), OUT_<name> (outgroups).
"""
import sys, re, collections
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from Bio import Phylo

tf, out = sys.argv[1], sys.argv[2]
mode = sys.argv[3] if len(sys.argv) > 3 else "collapsed"
INK, MUTED, BG = "#1b1b1b", "#6b6b6b", "#fdfdfb"
C = dict(adr1="#0F8F7E", nrg1="#0072B2", tir="#D17A00", ref="#C9184A", nic="#7A5195", out="#6b6b6b")
plt.rcParams.update({"font.family": "DejaVu Sans", "pdf.fonttype": 42, "svg.fonttype": "none"})

t = Phylo.read(tf, "newick")
leaves = t.get_terminals()
NRG1_L = {l.name for l in leaves if l.name.startswith("Callicarpa_nudiflora|CM111845.1-g274")}
TIR_L = {l.name for l in leaves if re.match(r"(Callicarpa_nudiflora\|CM111844|Elsholtzia_splendens\|CM112684)", l.name)}
def kind(n):
    if n in NRG1_L: return "nrg1"
    if n in TIR_L: return "tir"
    if n.startswith("OUT_"): return "out"
    if n.startswith("REF_"): return "nic" if n in ("REF_NbADR1", "REF_NRG1") else "ref"
    return "adr1"
adr1 = [l for l in leaves if kind(l.name) == "adr1"]
n_sp = len({l.name.split("|")[0] for l in adr1})
mrca = t.common_ancestor(adr1)
collapse = (mode == "collapsed")

# depths
depth = {t.root: 0.0}
for c in t.find_clades(order="preorder"):
    for ch in c.clades: depth[ch] = depth[c] + (ch.branch_length or 0)
xmax_tip = max(depth[l] for l in leaves)
pos = {}; cnt = [0]
def layout(c):
    if collapse and c is mrca:
        pos[c] = cnt[0] + 1.0; cnt[0] += 3; return pos[c]
    if c.is_terminal():
        pos[c] = cnt[0]; cnt[0] += 1; return pos[c]
    ys = [layout(ch) for ch in c.clades]; pos[c] = sum(ys) / len(ys); return pos[c]
layout(t.root)
H = cnt[0]

def sup(c):
    m = re.match(r"([\d.]+)/([\d.]+)", c.name or "")
    return (float(m.group(1)), float(m.group(2))) if m else None

fig_w, fig_h = (9.0, 6.4) if collapse else (9.5, 12.5)
fig, ax = plt.subplots(figsize=(fig_w, fig_h), facecolor=BG)
ax.set_facecolor(BG)
LW = 0.9
def draw(c):
    if collapse and c is mrca: return
    for ch in c.clades:
        ax.plot([depth[c], depth[ch]], [pos[ch]] * 2, c="#444", lw=LW, solid_capstyle="butt")
        draw(ch)
    if c.clades:
        ys = [pos[ch] for ch in c.clades]
        ax.plot([depth[c]] * 2, [min(ys), max(ys)], c="#444", lw=LW)
draw(t.root)
# stem to collapsed node
xlab = xmax_tip * 1.02
for l in leaves:
    if collapse and l in set(adr1): continue
    k = kind(l.name); nm = l.name
    if nm.startswith(("REF_", "OUT_")):
        nm = nm.split("_", 1)[1].replace("_", "-").replace("AtNRG1-", "AtNRG1.")
        if nm == "ADR1": nm = "ADR1†"
    else:
        sp, g = nm.split("|"); nm = sp.replace("_", " ") + "  " + g.split("-")[-1].replace(".t1", "")
    ax.text(depth[l] + xmax_tip * 0.012, pos[l], nm, va="center", ha="left", fontsize=7.5 if collapse else 6.4,
            color=C[k], style="italic" if k in ("adr1", "nrg1", "tir") else "normal", fontweight="bold" if k in ("ref", "nic", "out") else "normal")
if collapse:
    x0 = depth[mrca]; x1 = max(depth[l] for l in adr1); y = pos[mrca]
    ax.add_patch(Polygon([[x0, y], [x1, y - 1.25], [x1, y + 1.25]], closed=True, fc=C["adr1"], ec=C["adr1"], alpha=0.85, lw=0.8, zorder=3))
    ax.text(x1 + xmax_tip * 0.012, y, f"Lamiaceae ADR1\n{len(adr1)} proteins, {n_sp} species", va="center", fontsize=8.5, color=C["adr1"], fontweight="bold")
# support dots: strong = SH-aLRT>=90 and UFBoot>=95
for c in t.get_nonterminals():
    if collapse and c is not mrca:
        inside = set(adr1)
        if set(c.get_terminals()) <= inside: continue
    s = sup(c)
    if not s or c is t.root: continue
    strong = s[0] >= 90 and s[1] >= 95
    ax.scatter([depth[c]], [pos[c]], s=11 if collapse else 7, fc=INK if strong else BG, ec=INK, lw=0.7, zorder=5)
sb = max(0.1, round(xmax_tip / 8, 1))
ax.plot([0, sb], [-1.2, -1.2], c=INK, lw=1.2); ax.text(sb / 2, -1.9, f"{sb:g} substitutions/site", ha="center", fontsize=7, color=MUTED)
ax.set_xlim(-xmax_tip * 0.01, xmax_tip * (1.55 if collapse else 1.42)); ax.set_ylim(H + 0.5, -2.8); ax.axis("off")
ax.set_title("ADR1/NRG1 (RNL) lineage in Lamiaceae: maximum-likelihood tree of NB-ARC Clade 18", fontsize=10.5, color=INK, loc="left", pad=8)
from matplotlib.lines import Line2D
H_ = [Line2D([], [], marker="s", ls="", color=C["adr1"], label="Lamiaceae ADR1"), Line2D([], [], marker="s", ls="", color=C["nrg1"], label="Lamiaceae NRG1"),
      Line2D([], [], marker="s", ls="", color=C["tir"], label="Lamiaceae, probable TIR-NLR"), Line2D([], [], marker="s", ls="", color=C["ref"], label="reference, A. thaliana"),
      Line2D([], [], marker="s", ls="", color=C["nic"], label="reference, N. benthamiana"), Line2D([], [], marker="s", ls="", color=C["out"], label="outgroup (ZAR1, NRC0)")]
fig.legend(handles=H_, loc="lower center", ncol=3, frameon=False, fontsize=7.5, bbox_to_anchor=(0.5, 0.05))
fig.text(0.5, 0.012, "Q.PLANT+G4, IQ-TREE 3, 1,000 ultrafast bootstraps. Filled circle: SH-aLRT ≥ 90 and UFBoot ≥ 95; open circle: lower support. Tree drawn at the outgroup root.\n"
         "† ADR1 reference: host species not listed in the Table S1 copy used here (A. thaliana ADR1 expected).", ha="center", fontsize=6.4, color=MUTED)
fig.subplots_adjust(left=0.03, right=0.98, top=0.93, bottom=0.14 if collapse else 0.10)
fig.savefig(out + ".png", dpi=300, facecolor=BG); fig.savefig(out + ".pdf", facecolor=BG)
print("ok", mode, len(adr1), "ADR1 proteins,", n_sp, "species")
