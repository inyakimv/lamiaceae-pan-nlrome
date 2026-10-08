"""
Map HyPhy FEL/MEME positively-selected codon sites onto protein domain
boundaries (NB-ARC, LRR, TIR) for the three ancestral orthogroups shared
with Arabidopsis thaliana (OG0000003_arab, OG0000006_arab,
OG0000115_arab).

Rationale
---------
FEL and MEME report selection at the level of ALIGNMENT COLUMNS (codon
sites), which is consistent across all sequences in a given alignment
regardless of individual gaps. To interpret WHERE in the protein
selection is acting, we need to know which alignment columns fall
within which structural domain.

Strategy
--------
1. For each orthogroup, extract the Arabidopsis thaliana sequence from
   the codon alignment (it is present in all three groups and has an
   unambiguous, independently characterized domain architecture).
2. Translate it, strip gaps, and run hmmscan against the existing
   reduced Pfam profile set (NB-ARC, LRR, TIR) already prepared on the
   server, to get exact domain boundaries in UNGAPPED PROTEIN
   coordinates.
3. Build a lookup table: alignment column -> domain (or "other"), using
   the Arabidopsis sequence's own gap pattern to convert ungapped
   protein positions to alignment column indices.
4. Parse each orthogroup's FEL and MEME JSON output to get the list of
   significant sites (p <= 0.1) and their alignment column index.
5. Cross-reference: for each significant site, look up which domain
   (if any) that alignment column falls into, and report a summary
   table (site counts by domain, per orthogroup and per method).

Usage (run from ~/project on the server):
    conda activate ANNEVO   # for hmmscan availability, or use base env
                             # if HMMER is installed system-wide there too
    python3 map_selection_to_domains.py

Requirements: Biopython, HMMER (hmmscan) available in PATH, and the
existing reduced Pfam profile file at
    nlr_domain_search/nlr_profiles.hmm
"""

import json
import re
import subprocess
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_DIR = Path.home() / "project"
HYPHY_DIR = PROJECT_DIR / "05_selection_hyphy"
PROFILES_HMM = PROJECT_DIR / "nlr_domain_search" / "nlr_profiles.hmm"

ORTHOGROUPS = ["OG0000003_arab", "OG0000006_arab", "OG0000115_arab"]

# Domain profile name -> human-readable category
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
    """Return (record_id, aligned_nt_seq) for the Arabidopsis sequence
    in a codon alignment fasta."""
    for rec in SeqIO.parse(codon_alignment_path, "fasta"):
        if "Athaliana" in rec.id or "Arabidopsis" in rec.id:
            return rec.id, str(rec.seq)
    return None, None


def translate_ungapped(aligned_nt_seq):
    """Strip gaps from an aligned nucleotide sequence and translate to
    protein, returning (ungapped_nt, protein)."""
    ungapped_nt = aligned_nt_seq.replace("-", "").replace("N", "")
    # Trim to a multiple of 3 just in case
    trim = len(ungapped_nt) % 3
    if trim:
        ungapped_nt = ungapped_nt[:-trim]
    protein = str(Seq(ungapped_nt).translate(to_stop=True))
    return ungapped_nt, protein


def run_hmmscan_on_protein(protein_seq, out_domtbl):
    """Write the protein to a temp fasta and run hmmscan against the
    reduced Pfam profile set, producing a domtblout file."""
    tmp_fasta = out_domtbl.with_suffix(".fasta")
    with open(tmp_fasta, "w") as f:
        f.write(">query\n" + protein_seq + "\n")

    subprocess.run(
        ["hmmscan", "--domtblout", str(out_domtbl), str(PROFILES_HMM), str(tmp_fasta)],
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=True
    )


def parse_domtblout(domtbl_path):
    """Return a list of (domain_category, ali_from, ali_to) in
    UNGAPPED PROTEIN coordinates (1-based, inclusive), for domains
    recognized in DOMAIN_CATEGORY."""
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


def build_column_to_domain_map(aligned_nt_seq, domain_hits):
    """Build a dict: alignment CODON column (1-based) -> domain category
    (or None), using the ungapped protein position of each codon of the
    reference (Arabidopsis) sequence."""
    codon_col_to_domain = {}
    ungapped_codon_index = 0  # 1-based ungapped codon / amino acid position

    n_codons = len(aligned_nt_seq) // 3
    for col in range(1, n_codons + 1):
        codon = aligned_nt_seq[(col - 1) * 3: col * 3]
        if "-" in codon or codon == "":
            # Gap in the reference at this alignment column: cannot
            # assign a domain via the reference sequence.
            codon_col_to_domain[col] = None
            continue
        ungapped_codon_index += 1
        domain = None
        for cat, start, end in domain_hits:
            if start <= ungapped_codon_index <= end:
                domain = cat
                break
        codon_col_to_domain[col] = domain
    return codon_col_to_domain


def parse_fel_significant_sites(fel_json_path, pvalue_threshold=0.1):
    """Return list of (site_index_1based, p_value, direction) for FEL
    significant sites. direction is 'positive' or 'negative'."""
    with open(fel_json_path) as f:
        data = json.load(f)

    headers = [h[0] for h in data["MLE"]["headers"]]
    try:
        alpha_idx = headers.index("alpha")
        beta_idx = headers.index("beta")
        p_idx = headers.index("p-value")
    except ValueError:
        # Header names can vary slightly between HyPhy versions; fall back
        # to commonly known positions if exact names are not found.
        alpha_idx, beta_idx, p_idx = 0, 1, 4

    rows = data["MLE"]["content"]["0"]
    sig_sites = []
    for i, row in enumerate(rows, start=1):
        p = row[p_idx]
        if p is not None and p <= pvalue_threshold:
            direction = "positive" if row[beta_idx] > row[alpha_idx] else "negative"
            sig_sites.append((i, p, direction))
    return sig_sites


def parse_meme_significant_sites(meme_json_path, pvalue_threshold=0.1):
    """Return list of (site_index_1based, p_value) for MEME episodic
    diversifying selection significant sites."""
    with open(meme_json_path) as f:
        data = json.load(f)

    headers = [h[0] for h in data["MLE"]["headers"]]
    try:
        p_idx = headers.index("p-value")
    except ValueError:
        p_idx = 6  # common fallback position in MEME output

    rows = data["MLE"]["content"]["0"]
    sig_sites = []
    for i, row in enumerate(rows, start=1):
        p = row[p_idx]
        if p is not None and p <= pvalue_threshold:
            sig_sites.append((i, p))
    return sig_sites


def summarize_domain_hits(sig_sites, column_to_domain):
    """Given a list of (site_index, ...) and the column->domain map,
    return counts per domain category."""
    counts = {"NB-ARC": 0, "LRR": 0, "TIR": 0, "other/unmapped": 0}
    for entry in sig_sites:
        site = entry[0]
        domain = column_to_domain.get(site)
        if domain in counts:
            counts[domain] += 1
        else:
            counts["other/unmapped"] += 1
    return counts


def main():
    all_results = []

    for og in ORTHOGROUPS:
        print("=" * 70)
        print(f">> Processing {og}")
        print("=" * 70)

        og_dir = HYPHY_DIR / og
        codon_alignment_path = og_dir / f"{og}_codon_alignment.fasta"
        fel_json_path = og_dir / f"{og}_FEL.json"
        meme_json_path = og_dir / f"{og}_MEME.json"

        if not codon_alignment_path.exists():
            print(f"  [!] Missing codon alignment for {og}, skipping.")
            continue

        ref_id, ref_aligned_nt = get_arabidopsis_alignment_seq(codon_alignment_path)
        if ref_aligned_nt is None:
            print(f"  [!] No Arabidopsis sequence found in {og} alignment, skipping.")
            continue
        print(f"  Reference sequence: {ref_id}")

        ungapped_nt, protein = translate_ungapped(ref_aligned_nt)
        print(f"  Reference protein length (ungapped): {len(protein)} aa")

        domtbl_path = og_dir / f"{og}_reference_domains.domtblout"
        run_hmmscan_on_protein(protein, domtbl_path)
        domain_hits = parse_domtblout(domtbl_path)
        print(f"  Domain hits found: {domain_hits}")

        column_to_domain = build_column_to_domain_map(ref_aligned_nt, domain_hits)

        fel_sites = []
        if fel_json_path.exists():
            fel_sites = parse_fel_significant_sites(fel_json_path)
            fel_positive = [s for s in fel_sites if s[2] == "positive"]
            fel_counts = summarize_domain_hits(fel_positive, column_to_domain)
            print(f"  FEL positive sites (n={len(fel_positive)}): {fel_counts}")
        else:
            fel_counts = {}
            print("  [!] FEL JSON not found.")

        meme_counts = {}
        if meme_json_path.exists():
            meme_sites = parse_meme_significant_sites(meme_json_path)
            meme_counts = summarize_domain_hits(meme_sites, column_to_domain)
            print(f"  MEME episodic sites (n={len(meme_sites)}): {meme_counts}")
        else:
            print("  [!] MEME JSON not found.")

        all_results.append({
            "orthogroup": og,
            "reference_protein_length_aa": len(protein),
            "domain_hits": domain_hits,
            "fel_positive_by_domain": fel_counts,
            "meme_episodic_by_domain": meme_counts,
        })
        print()

    # --- Save summary as CSV ---
    import csv
    out_csv = PROJECT_DIR / "selection_by_domain_summary.csv"
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "orthogroup", "method", "NB-ARC", "LRR", "TIR", "other_unmapped"
        ])
        for r in all_results:
            writer.writerow([
                r["orthogroup"], "FEL_positive",
                r["fel_positive_by_domain"].get("NB-ARC", ""),
                r["fel_positive_by_domain"].get("LRR", ""),
                r["fel_positive_by_domain"].get("TIR", ""),
                r["fel_positive_by_domain"].get("other/unmapped", ""),
            ])
            writer.writerow([
                r["orthogroup"], "MEME_episodic",
                r["meme_episodic_by_domain"].get("NB-ARC", ""),
                r["meme_episodic_by_domain"].get("LRR", ""),
                r["meme_episodic_by_domain"].get("TIR", ""),
                r["meme_episodic_by_domain"].get("other/unmapped", ""),
            ])

    print("=" * 70)
    print(f"Summary saved to: {out_csv}")
    print("=" * 70)


if __name__ == "__main__":
    main()
