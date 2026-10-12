"""
Statistical analysis of Foldseek structural comparisons, for both the
factorial-design orthogroups (15 structures, recent lineage-specific
expansion) and the ancestral-core orthogroups (6 structures, anchored
by Arabidopsis thaliana orthologs).

No figures are generated here -- this script only computes numbers and
saves CSVs, to be used later when building the paper's figures.

Part 1 -- Factorial design (15 structures):
    Reads 08_structure_comparison/structure_comparison.tsv and computes,
    for every pairwise comparison, whether the pair belongs to the SAME
    orthogroup ("within-group") or to DIFFERENT orthogroups
    ("between-group"), based on the "ogXXXXXXX_" prefix in each
    structure's filename. Reports summary statistics for both groups
    and a Mann-Whitney U test comparing them (mirroring the approach
    already used for gene size distribution).

Part 2 -- Ancestral core (6 structures):
    Reads 08_structure_comparison_ancestral/structure_comparison_ancestral.tsv
    and extracts specifically the three Arabidopsis-vs-Lamiaceae pairs
    (one per orthogroup: OG0000003, OG0000006, OG0000115), reporting
    TM-score alongside sequence identity for each. No formal
    significance test is computed here (n=3 is too small), just a
    clean summary table.

Usage (run from ~/project):
    python3 foldseek_statistics.py

Requirements: pandas, scipy
"""

import re
import pandas as pd
from pathlib import Path
from scipy.stats import mannwhitneyu

PROJECT_DIR = Path.home() / "project"

# --- Part 1: factorial design ---
FACTORIAL_TSV = PROJECT_DIR / "08_structure_comparison" / "structure_comparison.tsv"
FACTORIAL_COLUMNS = ["query", "target", "fident", "alnlen", "evalue", "bits",
                     "qtmscore", "ttmscore", "alntmscore"]

# --- Part 2: ancestral core ---
ANCESTRAL_TSV = PROJECT_DIR / "08_structure_comparison_ancestral" / "structure_comparison_ancestral.tsv"


def extract_orthogroup(name):
    """Extract the ogXXXXXXX identifier from a structure filename such
    as 'og0000000_nc062966_1_nlr114' -> 'og0000000'."""
    m = re.match(r"(og\d+)_", name, re.IGNORECASE)
    return m.group(1).lower() if m else None


def describe(values):
    if not values:
        return {}
    s = pd.Series(values)
    return {
        "n": len(s),
        "mean": round(s.mean(), 3),
        "median": round(s.median(), 3),
        "stdev": round(s.std(), 3) if len(s) > 1 else 0,
        "min": round(s.min(), 3),
        "max": round(s.max(), 3),
    }


def analyze_factorial():
    if not FACTORIAL_TSV.exists():
        print(f"[!] {FACTORIAL_TSV} not found, skipping factorial design analysis.")
        return

    df = pd.read_csv(FACTORIAL_TSV, sep="\t", names=FACTORIAL_COLUMNS)

    # Remove self-comparisons (query == target, identity/TM-score == 1)
    df = df[df["query"] != df["target"]].copy()

    df["query_og"] = df["query"].apply(extract_orthogroup)
    df["target_og"] = df["target"].apply(extract_orthogroup)
    df = df.dropna(subset=["query_og", "target_og"])

    df["comparison_type"] = df.apply(
        lambda r: "within_group" if r["query_og"] == r["target_og"] else "between_group",
        axis=1
    )

    within = df[df["comparison_type"] == "within_group"]["alntmscore"].tolist()
    between = df[df["comparison_type"] == "between_group"]["alntmscore"].tolist()

    print("=" * 70)
    print("PART 1: FACTORIAL DESIGN (recent lineage-specific expansion)")
    print("=" * 70)
    print(f"Total pairwise comparisons (excluding self): {len(df)}")
    print(f"  Within-group comparisons: {len(within)}")
    print(f"  Between-group comparisons: {len(between)}")
    print()

    within_stats = describe(within)
    between_stats = describe(between)

    print(f"{'Metric':<10}{'Within-group':>15}{'Between-group':>17}")
    for key in ["n", "mean", "median", "stdev", "min", "max"]:
        print(f"{key:<10}{within_stats.get(key, ''):>15}{between_stats.get(key, ''):>17}")

    if len(within) > 0 and len(between) > 0:
        u_stat, p_value = mannwhitneyu(within, between, alternative="two-sided")
        print()
        print(f"Mann-Whitney U test (within vs. between): U={u_stat:.1f}, p={p_value:.4g}")

    # Save detailed per-pair CSV
    out_csv = PROJECT_DIR / "foldseek_factorial_within_between.csv"
    df[["query", "target", "query_og", "target_og", "fident",
        "alntmscore", "comparison_type"]].to_csv(out_csv, index=False)
    print(f"\nDetailed per-pair data saved to: {out_csv}")

    # Save summary CSV
    summary_rows = []
    for label, stats in [("within_group", within_stats), ("between_group", between_stats)]:
        row = {"comparison_type": label}
        row.update(stats)
        summary_rows.append(row)
    summary_out = PROJECT_DIR / "foldseek_factorial_summary.csv"
    pd.DataFrame(summary_rows).to_csv(summary_out, index=False)
    print(f"Summary statistics saved to: {summary_out}")
    print()


def analyze_ancestral():
    if not ANCESTRAL_TSV.exists():
        print(f"[!] {ANCESTRAL_TSV} not found, skipping ancestral core analysis.")
        return

    df = pd.read_csv(ANCESTRAL_TSV, sep="\t", names=FACTORIAL_COLUMNS)
    df = df[df["query"] != df["target"]].copy()

    # Identify Arabidopsis-vs-Lamiaceae pairs within the same orthogroup
    def og_of(name):
        m = re.match(r"fold_(og\d+)_", name, re.IGNORECASE)
        return m.group(1).lower() if m else None

    def is_athaliana(name):
        return "athaliana" in name.lower()

    df["query_og"] = df["query"].apply(og_of)
    df["target_og"] = df["target"].apply(og_of)
    df["query_is_ara"] = df["query"].apply(is_athaliana)
    df["target_is_ara"] = df["target"].apply(is_athaliana)

    pairs = df[
        (df["query_og"] == df["target_og"]) &
        (df["query_is_ara"] != df["target_is_ara"])
    ].copy()

    # Deduplicate (Foldseek reports both query->target and target->query)
    pairs["og_key"] = pairs["query_og"]
    pairs = pairs.drop_duplicates(subset=["og_key"], keep="first")

    print("=" * 70)
    print("PART 2: ANCESTRAL CORE (Arabidopsis-vs-Lamiaceae, per orthogroup)")
    print("=" * 70)
    print(f"{'Orthogroup':<12}{'Seq. identity':>15}{'alntmscore':>14}")
    for _, row in pairs.sort_values("og_key").iterrows():
        print(f"{row['og_key']:<12}{row['fident']:>15.3f}{row['alntmscore']:>14.3f}")

    print()
    print(f"Mean TM-score across the 3 ancestral pairs: {pairs['alntmscore'].mean():.3f}")
    print(f"Mean sequence identity across the 3 ancestral pairs: {pairs['fident'].mean():.3f}")
    print("(Note: n=3, too small for a formal significance test -- reported as")
    print(" descriptive summary only.)")

    out_csv = PROJECT_DIR / "foldseek_ancestral_summary.csv"
    pairs[["og_key", "query", "target", "fident", "alntmscore"]].rename(
        columns={"og_key": "orthogroup"}
    ).to_csv(out_csv, index=False)
    print(f"\nSummary saved to: {out_csv}")
    print()


if __name__ == "__main__":
    analyze_factorial()
    analyze_ancestral()
