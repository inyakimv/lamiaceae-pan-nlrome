"""
Extended analysis: list the exact protein position (Arabidopsis
reference coordinates) of every significant FEL/MEME site, for the
three ancestral orthogroups, so we can inspect whether selected sites
within a domain cluster in a specific subregion (e.g., near the
domain's distal/interface end) rather than being spread uniformly
across the whole domain (which would be more consistent with a
functional, interpretable pattern rather than noise).

This reuses the same domain-mapping logic as map_selection_to_domains.py
but outputs per-site detail instead of only aggregate counts.

Usage (run from ~/project, same environment as map_selection_to_domains.py):
    python3 list_selected_sites_detail.py
"""

import json
import csv
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_DIR = Path.home() / "project"
HYPHY_DIR = PROJECT_DIR / "05_selection_hyphy"

ORTHOGROUPS = ["OG0000003_arab", "OG0000006_arab", "OG0000115_arab"]

DOMAIN_CATEGORY = {
    "NB-ARC": "NB-ARC",
    "LRR_1": "LRR",
    "LRR_3": "LRR",
    "LRR_8": "LRR",
    "LRR_4": "LRR",
    "TIR": "TIR",
    "TIR_2": "TIR",
}


def get_arabidopsis_alignment_seq(codon_alignment_path):
    for rec in SeqIO.parse(codon_alignment_path, "fasta"):
        if "Athaliana" in rec.id or "Arabidopsis" in rec.id:
            return rec.id, str(rec.seq)
    return None, None


def parse_domtblout(domtbl_path):
    hits = []
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            profile_name = fields[0]
            if profile_name not in DOMAIN_CATEGORY:
                continue
            ali_from = int(fields[17])
            ali_to = int(fields[18])
            hits.append((DOMAIN_CATEGORY[profile_name], ali_from, ali_to))
    return hits


def build_column_maps(aligned_nt_seq, domain_hits):
    """Return two dicts: col -> arabidopsis_ungapped_protein_position
    (or None if reference has a gap there), and col -> domain (or None)."""
    col_to_pos = {}
    col_to_domain = {}
    ungapped_codon_index = 0
    n_codons = len(aligned_nt_seq) // 3

    for col in range(1, n_codons + 1):
        codon = aligned_nt_seq[(col - 1) * 3: col * 3]
        if "-" in codon or codon == "":
            col_to_pos[col] = None
            col_to_domain[col] = None
            continue
        ungapped_codon_index += 1
        col_to_pos[col] = ungapped_codon_index
        domain = None
        for cat, start, end in domain_hits:
            if start <= ungapped_codon_index <= end:
                domain = cat
                break
        col_to_domain[col] = domain
    return col_to_pos, col_to_domain


def parse_fel_sites(fel_json_path, pvalue_threshold=0.1):
    with open(fel_json_path) as f:
        data = json.load(f)
    headers = [h[0] for h in data["MLE"]["headers"]]
    alpha_idx = headers.index("alpha") if "alpha" in headers else 0
    beta_idx = headers.index("beta") if "beta" in headers else 1
    p_idx = headers.index("p-value") if "p-value" in headers else 4

    rows = data["MLE"]["content"]["0"]
    sites = []
    for i, row in enumerate(rows, start=1):
        p = row[p_idx]
        if p is not None and p <= pvalue_threshold:
            direction = "positive" if row[beta_idx] > row[alpha_idx] else "negative"
            sites.append((i, p, direction))
    return sites


def parse_meme_sites(meme_json_path, pvalue_threshold=0.1):
    with open(meme_json_path) as f:
        data = json.load(f)
    headers = [h[0] for h in data["MLE"]["headers"]]
    p_idx = headers.index("p-value") if "p-value" in headers else 6

    rows = data["MLE"]["content"]["0"]
    sites = []
    for i, row in enumerate(rows, start=1):
        p = row[p_idx]
        if p is not None and p <= pvalue_threshold:
            sites.append((i, p))
    return sites


def main():
    out_rows = []

    for og in ORTHOGROUPS:
        print("=" * 70)
        print(f">> {og}")
        print("=" * 70)

        og_dir = HYPHY_DIR / og
        codon_alignment_path = og_dir / f"{og}_codon_alignment.fasta"
        fel_json_path = og_dir / f"{og}_FEL.json"
        meme_json_path = og_dir / f"{og}_MEME.json"
        domtbl_path = og_dir / f"{og}_reference_domains.domtblout"

        if not domtbl_path.exists():
            print(f"  [!] {domtbl_path} not found — run map_selection_to_domains.py first.")
            continue

        ref_id, ref_aligned_nt = get_arabidopsis_alignment_seq(codon_alignment_path)
        domain_hits = parse_domtblout(domtbl_path)
        col_to_pos, col_to_domain = build_column_maps(ref_aligned_nt, domain_hits)

        # FEL positive sites, with detail
        if fel_json_path.exists():
            fel_sites = parse_fel_sites(fel_json_path)
            fel_positive = [s for s in fel_sites if s[2] == "positive"]
            print(f"\n  FEL positive sites (n={len(fel_positive)}):")
            print(f"  {'AlnCol':<8}{'AraPos':<8}{'Domain':<10}{'p-value':<10}")
            for col, p, _ in fel_positive:
                pos = col_to_pos.get(col)
                dom = col_to_domain.get(col) or "-"
                print(f"  {col:<8}{str(pos):<8}{dom:<10}{p:<10.4f}")
                out_rows.append([og, "FEL_positive", col, pos, dom, round(p, 4)])

        # MEME episodic sites, with detail
        if meme_json_path.exists():
            meme_sites = parse_meme_sites(meme_json_path)
            print(f"\n  MEME episodic sites (n={len(meme_sites)}):")
            print(f"  {'AlnCol':<8}{'AraPos':<8}{'Domain':<10}{'p-value':<10}")
            for col, p in meme_sites:
                pos = col_to_pos.get(col)
                dom = col_to_domain.get(col) or "-"
                print(f"  {col:<8}{str(pos):<8}{dom:<10}{p:<10.4f}")
                out_rows.append([og, "MEME_episodic", col, pos, dom, round(p, 4)])

        print()

    out_csv = PROJECT_DIR / "selection_sites_detail.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["orthogroup", "method", "alignment_column",
                          "arabidopsis_protein_position", "domain", "p_value"])
        writer.writerows(out_rows)

    print("=" * 70)
    print(f"Full per-site detail saved to: {out_csv}")
    print("=" * 70)


if __name__ == "__main__":
    main()
