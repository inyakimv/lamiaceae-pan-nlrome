#!/usr/bin/env python3
"""
Reads the FEL and MEME JSON output files for one or more
orthogroups already processed by HyPhy, and produces two
summary files:

  1. hyphy_summary.csv        -- one row per orthogroup, with
                                  site counts (matches the
                                  numbers HyPhy prints in its
                                  own log files, used here as
                                  an internal sanity check)
  2. hyphy_significant_sites.csv -- one row per significant
                                  site, tagging whether it was
                                  found by FEL, MEME, or BOTH
                                  (the most robust candidates)

Usage:
    cd ~/project
    python3 summarize_hyphy_results_en.py

Automatically finds all */_FEL.json and */_MEME.json pairs
inside 05_selection_hyphy/.
"""

import json
import os
import glob
import csv

PROJECT_DIR = os.path.expanduser("~/project")
HYPHY_DIR = os.path.join(PROJECT_DIR, "05_selection_hyphy")
P_THRESHOLD = 0.1  # matches the threshold HyPhy itself used (p <= 0.1)


def find_header_index(headers, match_in="name", contains=None, exact=None):
    """headers: list of [name, description] pairs from HyPhy JSON.
    Finds the index of the column matching by name or description."""
    for i, h in enumerate(headers):
        name, desc = h[0], h[1] if len(h) > 1 else ""
        target = name if match_in == "name" else desc
        if exact is not None and target == exact:
            return i
        if contains is not None and contains.lower() in target.lower():
            return i
    return None


def parse_fel(json_path):
    """Returns list of (site_1based, p_value, direction)."""
    with open(json_path) as f:
        data = json.load(f)

    headers = data["MLE"]["headers"]
    rows = data["MLE"]["content"]["0"]

    idx_alpha = find_header_index(headers, "name", exact="alpha")
    idx_beta = find_header_index(headers, "name", exact="beta")
    idx_p = find_header_index(headers, "description", contains="p-value")

    if idx_p is None:
        raise ValueError(f"Could not find p-value column in {json_path}. Headers: {headers}")

    results = []
    for site_idx, row in enumerate(rows, start=1):
        p = row[idx_p]
        if p is None or p > P_THRESHOLD:
            continue
        direction = "unknown"
        if idx_alpha is not None and idx_beta is not None:
            direction = "positive" if row[idx_beta] > row[idx_alpha] else "negative"
        results.append((site_idx, p, direction))
    return results


def parse_meme(json_path):
    """Returns list of (site_1based, p_value). MEME only reports
    episodic POSITIVE selection, so no direction needed."""
    with open(json_path) as f:
        data = json.load(f)

    headers = data["MLE"]["headers"]
    rows = data["MLE"]["content"]["0"]

    idx_p = find_header_index(headers, "description", contains="p-value")
    if idx_p is None:
        raise ValueError(f"Could not find p-value column in {json_path}. Headers: {headers}")

    results = []
    for site_idx, row in enumerate(rows, start=1):
        p = row[idx_p]
        if p is None or p > P_THRESHOLD:
            continue
        results.append((site_idx, p))
    return results


def main():
    fel_files = sorted(glob.glob(os.path.join(HYPHY_DIR, "*", "*_FEL.json")))
    if not fel_files:
        print(f"ERROR: no *_FEL.json files found inside {HYPHY_DIR}")
        return

    summary_rows = []
    site_rows = []

    for fel_path in fel_files:
        og_dir = os.path.dirname(fel_path)
        og_id = os.path.basename(og_dir)
        meme_path = os.path.join(og_dir, f"{og_id}_MEME.json")

        print(f"Processing {og_id}...")

        fel_sites = parse_fel(fel_path)
        fel_positive = {s for s, p, d in fel_sites if d == "positive"}
        fel_negative = {s for s, p, d in fel_sites if d == "negative"}

        meme_sites = []
        meme_positive = set()
        if os.path.exists(meme_path):
            meme_sites = parse_meme(meme_path)
            meme_positive = {s for s, p in meme_sites}
        else:
            print(f"  [!] No MEME file found for {og_id}, skipping MEME columns.")

        overlap = fel_positive & meme_positive

        summary_rows.append({
            "orthogroup": og_id,
            "FEL_positive_sites": len(fel_positive),
            "FEL_negative_sites": len(fel_negative),
            "MEME_episodic_sites": len(meme_positive),
            "sites_significant_in_BOTH_FEL_and_MEME": len(overlap),
        })

        fel_p = {s: p for s, p, d in fel_sites}
        fel_dir = {s: d for s, p, d in fel_sites}
        meme_p = {s: p for s, p in meme_sites}

        all_sites = set(fel_p) | set(meme_p)
        for site in sorted(all_sites):
            in_fel = site in fel_p
            in_meme = site in meme_p
            method = "BOTH" if (in_fel and in_meme) else ("FEL" if in_fel else "MEME")
            site_rows.append({
                "orthogroup": og_id,
                "site": site,
                "method": method,
                "FEL_p_value": fel_p.get(site, ""),
                "FEL_direction": fel_dir.get(site, ""),
                "MEME_p_value": meme_p.get(site, ""),
            })

        print(f"  -> FEL: {len(fel_positive)} positive, {len(fel_negative)} negative | "
              f"MEME: {len(meme_positive)} episodic | Overlap (FEL+ & MEME): {len(overlap)}")

    summary_path = os.path.join(HYPHY_DIR, "hyphy_summary.csv")
    with open(summary_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "orthogroup", "FEL_positive_sites", "FEL_negative_sites",
            "MEME_episodic_sites", "sites_significant_in_BOTH_FEL_and_MEME"
        ])
        writer.writeheader()
        writer.writerows(summary_rows)

    sites_path = os.path.join(HYPHY_DIR, "hyphy_significant_sites.csv")
    with open(sites_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "orthogroup", "site", "method", "FEL_p_value", "FEL_direction", "MEME_p_value"
        ])
        writer.writeheader()
        writer.writerows(site_rows)

    print(f"\nSummary saved to: {summary_path}")
    print(f"Detailed significant sites saved to: {sites_path}")
    print(f"\nNOTE: compare the FEL/MEME counts above against each ")
    print(f"*_FEL.log / *_MEME.log 'Found ...' line as a sanity check.")


if __name__ == "__main__":
    main()
