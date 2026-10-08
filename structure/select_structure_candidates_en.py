#!/usr/bin/env python3
"""
Selects a small, prioritized set of candidate proteins to fold
with AlphaFold Server, based on the HyPhy selection results
already computed.

Strategy: for each of the 5 highlighted orthogroups, count how
many "BOTH" (FEL+MEME) significant sites fall within each
individual sequence's aligned region, and pick the top 2-3
sequences per group with the most "BOTH" sites -- these are
the most robust candidates for structural analysis.

Usage:
    cd ~/project
    python3 select_structure_candidates_en.py

Requires:
    - 05_selection_hyphy/hyphy_significant_sites.csv
    - 05_selection_hyphy/<OG>/<OG>_codon_alignment.fasta (to map
      alignment columns back to individual sequences)

Output:
    - 06_structure_retrieval/candidates_for_alphafold.fasta
      (ungapped protein sequences, ready to paste into
      AlphaFold Server)
    - 06_structure_retrieval/candidates_summary.csv
"""
import os
import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

PROJECT_DIR = os.path.expanduser("~/project")
HYPHY_DIR = os.path.join(PROJECT_DIR, "05_selection_hyphy")
OUT_DIR = os.path.join(PROJECT_DIR, "06_structure_retrieval")
N_PER_GROUP = 3  # how many top sequences to pick per orthogroup


def load_both_sites(sites_csv):
    """Returns dict: orthogroup -> set of alignment-column sites
    significant in BOTH FEL and MEME."""
    both = defaultdict(set)
    with open(sites_csv) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if row["method"] == "BOTH":
                both[row["orthogroup"]].add(int(row["site"]))
    return both


def count_covered_sites(aligned_seq, both_sites):
    """Counts how many of the BOTH-significant alignment columns
    are non-gap (i.e. actually present) in this sequence."""
    count = 0
    for site in both_sites:
        # site is 1-based column index in the alignment
        if site - 1 < len(aligned_seq) and aligned_seq[site - 1] != "-":
            count += 1
    return count


def main():
    sites_csv = os.path.join(HYPHY_DIR, "hyphy_significant_sites.csv")
    if not os.path.exists(sites_csv):
        print(f"ERROR: cannot find {sites_csv}")
        return

    both_by_og = load_both_sites(sites_csv)
    os.makedirs(OUT_DIR, exist_ok=True)

    all_candidates = []

    for og_id, both_sites in both_by_og.items():
        alignment_path = os.path.join(HYPHY_DIR, og_id, f"{og_id}_codon_alignment.fasta")
        if not os.path.exists(alignment_path):
            print(f"[!] No codon alignment found for {og_id}, skipping.")
            continue

        scored = []
        for rec in SeqIO.parse(alignment_path, "fasta"):
            # translate this codon-aligned sequence back to protein
            # (gaps '---' stay as '-' after translation for counting)
            nt = str(rec.seq)
            protein_chars = []
            for i in range(0, len(nt) - len(nt) % 3, 3):
                codon = nt[i:i+3]
                if codon == "---":
                    protein_chars.append("-")
                else:
                    protein_chars.append(str(Seq(codon).translate()))
            aligned_protein = "".join(protein_chars)

            n_covered = count_covered_sites(aligned_protein, both_sites)
            scored.append((n_covered, rec.id, aligned_protein.replace("-", "")))

        scored.sort(key=lambda x: -x[0])
        top = scored[:N_PER_GROUP]

        for n_covered, seq_id, ungapped_protein in top:
            all_candidates.append({
                "orthogroup": og_id,
                "sequence_id": seq_id,
                "n_both_sites_covered": n_covered,
                "protein_length": len(ungapped_protein),
                "protein_seq": ungapped_protein,
            })
            print(f"{og_id}: {seq_id} -> {n_covered} BOTH-sites covered, "
                  f"{len(ungapped_protein)} aa")

    # Write FASTA ready for AlphaFold Server
    fasta_out = os.path.join(OUT_DIR, "candidates_for_alphafold.fasta")
    with open(fasta_out, "w") as f:
        for c in all_candidates:
            short_name = f"{c['orthogroup']}_{c['sequence_id'].split('|')[-2] if '|' in c['sequence_id'] else c['sequence_id']}"
            f.write(f">{short_name}\n{c['protein_seq']}\n")

    # Write summary CSV
    csv_out = os.path.join(OUT_DIR, "candidates_summary.csv")
    with open(csv_out, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "orthogroup", "sequence_id", "n_both_sites_covered", "protein_length"
        ])
        writer.writeheader()
        for c in all_candidates:
            writer.writerow({k: c[k] for k in writer.fieldnames})

    print(f"\n{len(all_candidates)} candidates selected.")
    print(f"FASTA for AlphaFold Server: {fasta_out}")
    print(f"Summary: {csv_out}")


if __name__ == "__main__":
    main()
