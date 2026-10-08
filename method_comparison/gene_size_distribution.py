"""
Gene size distribution statistics: NLR-Annotator vs ANNEVO+HMMER
complete-like NLR loci, across the 46 Lamiaceae species.

Publication-ready version:
  - Mann-Whitney U test (non-parametric, appropriate for skewed data)
  - Colorblind-safe palette (Okabe-Ito)
  - Vector export (PDF) in addition to PNG
  - Caption text saved separately (title removed from the figure itself)
  - Note on non-independence of observations (genes from 46 species
    pooled together; phylogenetic structure not accounted for)

Usage:
    cd ~/project
    python3 gene_size_distribution.py

Requirements: pandas, matplotlib, scipy
    pip install scipy --break-system-packages
"""

import re
from pathlib import Path
from collections import defaultdict
import pandas as pd
import statistics as stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu

NLR_ANNOTATOR_DIR = Path("02_nlr_annotation")
ANNEVO_GFF_DIR = Path("00_genomes")
DOMTBL_DIR = Path("nlr_domain_search/per_species")

LRR_PROFILES = {"LRR_1", "LRR_3", "LRR_8", "LRR_4"}
COMPLETE_ARCHITECTURES = {"CC-NBARC-LRR", "TIR-NBARC-LRR"}


def load_nlr_annotator_lengths(genome):
    species_dir = NLR_ANNOTATOR_DIR / genome
    if not species_dir.is_dir():
        return None
    nlr_files = list(species_dir.glob("*_nlr.txt"))
    if not nlr_files:
        return None
    lengths = []
    with open(nlr_files[0]) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            arch, start, end = fields[2], fields[3], fields[4]
            if arch in COMPLETE_ARCHITECTURES:
                lengths.append(int(end) - int(start) + 1)
    return lengths


def load_annevo_complete_like_lengths(genome):
    gff_path = ANNEVO_GFF_DIR / f"{genome}.annevo.gff3"
    domtbl_path = DOMTBL_DIR / f"{genome}.domtblout"
    if not gff_path.exists() or not domtbl_path.exists():
        return None

    # Gene coordinates
    gene_coords = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            start, end, attrs = fields[3], fields[4], fields[8]
            m = re.search(r"ID=([^;]+)", attrs)
            if not m:
                continue
            gene_coords[m.group(1)] = (int(start), int(end))

    # Complete-like protein -> gene IDs
    protein_domains = defaultdict(set)
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            protein_domains[fields[3]].add(fields[0])

    lengths = []
    for protein_id, domains in protein_domains.items():
        if "NB-ARC" in domains and (domains & LRR_PROFILES):
            gene_id = re.sub(r"\.t\d+$", "", protein_id)
            if gene_id in gene_coords:
                start, end = gene_coords[gene_id]
                lengths.append(end - start + 1)
    return lengths


def describe(lengths):
    if not lengths:
        return {}
    sorted_l = sorted(lengths)
    n = len(sorted_l)
    return {
        "n": n,
        "mean": round(stats.mean(sorted_l), 1),
        "median": round(stats.median(sorted_l), 1),
        "stdev": round(stats.stdev(sorted_l), 1) if n > 1 else 0,
        "min": sorted_l[0],
        "max": sorted_l[-1],
        "q1": round(sorted_l[int(n * 0.25)], 1),
        "q3": round(sorted_l[int(n * 0.75)], 1),
    }


genomes = sorted(p.name for p in NLR_ANNOTATOR_DIR.iterdir() if p.is_dir())
print(f"Found {len(genomes)} species folders.\n")

all_nlr_annotator_lengths = []
all_annevo_lengths = []
per_species_rows = []

for genome in genomes:
    nlr_lengths = load_nlr_annotator_lengths(genome)
    annevo_lengths = load_annevo_complete_like_lengths(genome)

    if nlr_lengths is None or annevo_lengths is None:
        print(f"  [!] Skipping {genome}: missing input file(s).")
        continue

    all_nlr_annotator_lengths.extend(nlr_lengths)
    all_annevo_lengths.extend(annevo_lengths)

    nlr_stats = describe(nlr_lengths)
    annevo_stats = describe(annevo_lengths)

    per_species_rows.append({
        "genome": genome,
        "nlr_annotator_n": nlr_stats.get("n"),
        "nlr_annotator_mean_bp": nlr_stats.get("mean"),
        "nlr_annotator_median_bp": nlr_stats.get("median"),
        "annevo_n": annevo_stats.get("n"),
        "annevo_mean_bp": annevo_stats.get("mean"),
        "annevo_median_bp": annevo_stats.get("median"),
    })

# --- Save per-species summary ---
per_species_df = pd.DataFrame(per_species_rows).sort_values("genome")
per_species_df.to_csv("gene_size_summary_per_species.csv", index=False)

# --- Save raw lengths for plotting ---
pd.DataFrame({"length_bp": all_nlr_annotator_lengths, "method": "NLR-Annotator"}).to_csv(
    "gene_lengths_nlr_annotator.csv", index=False)
pd.DataFrame({"length_bp": all_annevo_lengths, "method": "ANNEVO+HMMER"}).to_csv(
    "gene_lengths_annevo.csv", index=False)

# --- Aggregate statistics ---
agg_nlr = describe(all_nlr_annotator_lengths)
agg_annevo = describe(all_annevo_lengths)

print("=" * 70)
print("AGGREGATE GENE SIZE STATISTICS (all 46 species combined)")
print("=" * 70)
print(f"{'Metric':<12}{'NLR-Annotator':>18}{'ANNEVO+HMMER':>18}")
for key in ["n", "mean", "median", "stdev", "min", "q1", "q3", "max"]:
    print(f"{key:<12}{agg_nlr.get(key, ''):>18}{agg_annevo.get(key, ''):>18}")

print()
print("Saved:")
print("  gene_size_summary_per_species.csv (per-species stats)")
print("  gene_lengths_nlr_annotator.csv (raw lengths, NLR-Annotator)")
print("  gene_lengths_annevo.csv (raw lengths, ANNEVO+HMMER)")
print("=" * 70)

# --- Statistical test: Mann-Whitney U ---
# Chosen over a t-test because the data are strongly right-skewed
# (long tail of large genes), violating the normality assumption of
# a t-test. Mann-Whitney U tests whether one distribution is
# stochastically greater than the other, without assuming normality.
u_stat, p_value = mannwhitneyu(all_nlr_annotator_lengths, all_annevo_lengths,
                                alternative="two-sided")

print()
print("MANN-WHITNEY U TEST (NLR-Annotator vs ANNEVO+HMMER gene sizes)")
print(f"  U statistic: {u_stat:.1f}")
print(f"  p-value: {p_value:.3e}")
print(f"  Significant at alpha=0.05: {'YES' if p_value < 0.05 else 'NO'}")
print("=" * 70)

# --- Important caveat, printed for the record ---
print()
print("CAVEAT (include in Methods/Discussion when reporting this test):")
print("  Observations pooled across 46 species are NOT fully independent")
print("  (genes from the same species, and species sharing phylogenetic")
print("  history, are not independent draws). The Mann-Whitney U result")
print("  above should be interpreted as descriptive evidence of a")
print("  systematic difference in gene size between methods, not as a")
print("  formally independent-samples test in the strict statistical")
print("  sense. A mixed-effects or per-species paired approach would be")
print("  needed for a fully rigorous test if this is required by reviewers.")
print("=" * 70)

# --- Generate violin plot comparing both distributions (log scale) ---
# Okabe-Ito colorblind-safe palette: blue (#0072B2) and orange (#E69F00)
fig, ax = plt.subplots(figsize=(8, 7))

violin_data = [all_nlr_annotator_lengths, all_annevo_lengths]
parts = ax.violinplot(violin_data, positions=[1, 2], showmedians=True, showextrema=True, widths=0.7)

colors = ["#0072B2", "#E69F00"]
for pc, color in zip(parts["bodies"], colors):
    pc.set_facecolor(color)
    pc.set_alpha(0.6)
    pc.set_edgecolor("black")
    pc.set_linewidth(0.8)

for part_name in ("cmedians", "cmins", "cmaxes", "cbars"):
    if part_name in parts:
        parts[part_name].set_edgecolor("black")
        parts[part_name].set_linewidth(1)

ax.set_yscale("log")
ax.set_xticks([1, 2])
ax.set_xticklabels([f"NLR-Annotator\n(n={len(all_nlr_annotator_lengths)})",
                     f"ANNEVO+HMMER\n(n={len(all_annevo_lengths)})"])
ax.set_ylabel("NLR gene size (bp, log scale)")
# NOTE: no in-figure title -- for publication, the descriptive title
# belongs in the figure caption (see gene_size_distribution_caption.txt),
# not baked into the image itself.
ax.yaxis.grid(True, which="major", linestyle="--", alpha=0.4)
ax.yaxis.grid(True, which="minor", linestyle=":", alpha=0.2)

# Annotate actual medians as text (offset to avoid overlapping the median line)
ax.text(1.28, agg_nlr["median"], f'median = {agg_nlr["median"]:.0f} bp',
        va="center", ha="left", fontsize=9, color="#004d80",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
ax.text(2.28, agg_annevo["median"], f'median = {agg_annevo["median"]:.0f} bp',
        va="center", ha="left", fontsize=9, color="#8a5a00",
        bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))

# Annotate Mann-Whitney U p-value with a significance bracket
y_max = max(max(all_nlr_annotator_lengths), max(all_annevo_lengths))
bracket_y = y_max * 1.6
ax.plot([1, 1, 2, 2], [bracket_y, bracket_y * 1.08, bracket_y * 1.08, bracket_y],
        color="black", linewidth=1)
p_label = "p < 0.001" if p_value < 0.001 else f"p = {p_value:.3f}"
ax.text(1.5, bracket_y * 1.12, p_label, ha="center", va="bottom", fontsize=9)

ax.set_xlim(0.5, 3.0)
ax.set_ylim(top=bracket_y * 2.2)

plt.tight_layout()
plt.savefig("gene_size_distribution.png", dpi=300)
plt.savefig("gene_size_distribution.pdf")
print("  gene_size_distribution.png (300 dpi, raster)")
print("  gene_size_distribution.pdf (vector, publication-ready)")

# --- Save figure caption separately (title lives here, not in the image) ---
caption = (
    "Figure X. NLR gene size distribution predicted by NLR-Annotator "
    "(motif-based genome scan) versus ANNEVO+HMMER (full gene annotation "
    "plus NB-ARC/LRR domain search), pooled across 46 Lamiaceae species. "
    f"NLR-Annotator: n={len(all_nlr_annotator_lengths)}, median="
    f"{agg_nlr['median']:.0f} bp. ANNEVO+HMMER: n={len(all_annevo_lengths)}, "
    f"median={agg_annevo['median']:.0f} bp. Distributions compared with a "
    f"two-sided Mann-Whitney U test (U={u_stat:.1f}, {p_label}). Y-axis on "
    "a log10 scale. Note: observations are pooled across species and are "
    "not fully statistically independent (see main text for discussion)."
)
with open("gene_size_distribution_caption.txt", "w") as f:
    f.write(caption)
print("  gene_size_distribution_caption.txt (figure caption)")
print("=" * 70)
