"""
Refined comparison between NLR-Annotator and ANNEVO+HMMER results.

Unlike the first pass (which only counted proteins with NB-ARC),
this script parses each per-species .domtblout file directly and
requires NB-ARC + at least one LRR domain in the SAME protein to
count it as "complete-like NLR" -- a fairer comparison against
NLR-Annotator's complete_NLR metric (which requires NB-ARC + LRR,
plus optionally TIR/CC).

It also flags TIR presence as a proxy for TIR-NBARC-LRR vs
CC-NBARC-LRR architecture (Pfam does not have a reliable CC-domain
profile comparable to NLR-Annotator's custom CC motif, so CC itself
cannot be scored directly here -- this is documented as a limitation).

Usage:
    cd ~/project
    python3 compare_nlr_methods_refined.py
"""

import pandas as pd
from pathlib import Path
from collections import defaultdict

DOMTBL_DIR = Path("nlr_domain_search/per_species")
NLR_ANNOTATOR_CSV = Path("02_nlr_annotation/nlr_summary.csv")
OUT_CSV = Path("nlr_method_comparison_detailed.csv")

LRR_PROFILES = {"LRR_1", "LRR_3", "LRR_8", "LRR_4"}
TIR_PROFILES = {"TIR", "TIR_2"}


def parse_domtblout(path):
    """Return, per protein ID, the set of HMM profile names that hit it."""
    protein_domains = defaultdict(set)
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            hmm_name = fields[0]      # target = HMM profile name
            protein_id = fields[3]    # query = protein ID
            protein_domains[protein_id].add(hmm_name)
    return protein_domains


rows = []
domtbl_files = sorted(DOMTBL_DIR.glob("*.domtblout"))
if not domtbl_files:
    raise SystemExit(f"ERROR: no .domtblout files found in {DOMTBL_DIR}")

print(f"Found {len(domtbl_files)} .domtblout files.\n")

for domtbl in domtbl_files:
    genome = domtbl.stem
    protein_domains = parse_domtblout(domtbl)

    n_nbarc = 0
    n_nbarc_lrr = 0          # NB-ARC + LRR in the same protein ("complete-like")
    n_nbarc_lrr_tir = 0      # NB-ARC + LRR + TIR ("TIR-NBARC-LRR-like")
    n_nbarc_lrr_no_tir = 0   # NB-ARC + LRR, no TIR ("CC/other-NBARC-LRR-like")

    for protein_id, domains in protein_domains.items():
        has_nbarc = "NB-ARC" in domains
        has_lrr = bool(domains & LRR_PROFILES)
        has_tir = bool(domains & TIR_PROFILES)

        if has_nbarc:
            n_nbarc += 1
            if has_lrr:
                n_nbarc_lrr += 1
                if has_tir:
                    n_nbarc_lrr_tir += 1
                else:
                    n_nbarc_lrr_no_tir += 1

    rows.append({
        "genome": genome,
        "annevo_hmmer_nbarc_only": n_nbarc,
        "annevo_hmmer_complete_like": n_nbarc_lrr,
        "annevo_hmmer_tir_nbarc_lrr_like": n_nbarc_lrr_tir,
        "annevo_hmmer_other_nbarc_lrr_like": n_nbarc_lrr_no_tir,
    })

hmmer_df = pd.DataFrame(rows).sort_values("genome")

# --- Load NLR-Annotator results ---
nlr_annotator = pd.read_csv(NLR_ANNOTATOR_CSV)
nlr_annotator = nlr_annotator.rename(columns={
    "genoma": "genome",
    "total_candidatos": "nlr_annotator_total",
    "CC_NBARC_LRR": "nlr_annotator_cc_nbarc_lrr",
    "TIR_NBARC_LRR": "nlr_annotator_tir_nbarc_lrr",
    "NLR_completos": "nlr_annotator_complete",
})[[
    "genome", "nlr_annotator_total", "nlr_annotator_cc_nbarc_lrr",
    "nlr_annotator_tir_nbarc_lrr", "nlr_annotator_complete"
]]

# --- Merge ---
merged = pd.merge(nlr_annotator, hmmer_df, on="genome", how="outer", indicator=True)

mismatches = merged[merged["_merge"] != "both"]
if len(mismatches) > 0:
    print("WARNING: genome name mismatches between the two datasets:")
    print(mismatches[["genome", "_merge"]].to_string(index=False))
    print()

merged = merged.drop(columns=["_merge"])

# --- Comparable metric: NLR-Annotator complete vs HMMER complete-like ---
merged["difference"] = merged["annevo_hmmer_complete_like"] - merged["nlr_annotator_complete"]
merged["ratio_hmmer_vs_annotator"] = (
    merged["annevo_hmmer_complete_like"] / merged["nlr_annotator_complete"]
).round(2)

merged = merged.sort_values("genome")
merged.to_csv(OUT_CSV, index=False)

# --- Summary ---
print("=" * 60)
print(f"Total genomes compared: {len(merged)}")
print()
print("Totals across all 46 species:")
print(f"  NLR-Annotator complete NLRs:              {merged['nlr_annotator_complete'].sum()}")
print(f"  ANNEVO+HMMER NB-ARC only (any):            {merged['annevo_hmmer_nbarc_only'].sum()}")
print(f"  ANNEVO+HMMER complete-like (NB-ARC+LRR):   {merged['annevo_hmmer_complete_like'].sum()}")
print(f"  ANNEVO+HMMER TIR-NBARC-LRR-like:           {merged['annevo_hmmer_tir_nbarc_lrr_like'].sum()}")
print(f"  ANNEVO+HMMER other-NBARC-LRR-like:         {merged['annevo_hmmer_other_nbarc_lrr_like'].sum()}")
print()
print("Comparison (NLR-Annotator complete vs ANNEVO+HMMER complete-like):")
print(f"  Mean ratio (HMMER/NLR-Annotator):  {merged['ratio_hmmer_vs_annotator'].mean():.2f}")
print(f"  Correlation (Pearson r):           {merged['nlr_annotator_complete'].corr(merged['annevo_hmmer_complete_like']):.3f}")
print()
print(f"Saved detailed per-species results to: {OUT_CSV}")
print("=" * 60)
