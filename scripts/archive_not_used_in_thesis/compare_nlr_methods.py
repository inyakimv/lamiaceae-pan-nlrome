import pandas as pd

nlr_annotator = pd.read_csv("02_nlr_annotation/nlr_summary.csv")
nlr_annotator = nlr_annotator.rename(columns={
    "genoma": "genome",
    "total_candidatos": "nlr_annotator_total",
    "NLR_completos": "nlr_annotator_complete"
})[["genome", "nlr_annotator_total", "nlr_annotator_complete"]]

hmmer = pd.read_csv("nlr_domain_search/hmmer_nlr_summary.csv")
hmmer = hmmer.rename(columns={
    "proteins_with_NBARC": "annevo_hmmer_nbarc_proteins"
})[["genome", "annevo_hmmer_nbarc_proteins"]]

merged = pd.merge(nlr_annotator, hmmer, on="genome", how="outer", indicator=True)

mismatches = merged[merged["_merge"] != "both"]
if len(mismatches) > 0:
    print("WARNING: genome name mismatches between the two datasets:")
    print(mismatches[["genome", "_merge"]].to_string(index=False))
    print()

merged = merged.drop(columns=["_merge"])
merged["difference"] = merged["annevo_hmmer_nbarc_proteins"] - merged["nlr_annotator_complete"]
merged["ratio_hmmer_vs_annotator"] = (merged["annevo_hmmer_nbarc_proteins"] / merged["nlr_annotator_complete"]).round(2)
merged = merged.sort_values("genome")
merged.to_csv("nlr_method_comparison.csv", index=False)

print(f"Total genomes compared: {len(merged)}")
print()
print("Summary statistics:")
print(f"  NLR-Annotator complete NLRs (total): {merged['nlr_annotator_complete'].sum()}")
print(f"  ANNEVO+HMMER NB-ARC proteins (total): {merged['annevo_hmmer_nbarc_proteins'].sum()}")
print(f"  Mean ratio (HMMER/NLR-Annotator): {merged['ratio_hmmer_vs_annotator'].mean():.2f}")
print(f"  Correlation (Pearson r): {merged['nlr_annotator_complete'].corr(merged['annevo_hmmer_nbarc_proteins']):.3f}")
print()
print("Saved to: nlr_method_comparison.csv")
