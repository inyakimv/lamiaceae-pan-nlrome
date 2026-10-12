"""
Locus-level (genomic coordinate) comparison between NLR-Annotator and
ANNEVO+HMMER complete-like NLR candidates.

For each of the 46 species, this script:
  1. Reads NLR-Annotator's complete candidates (architecture ==
     CC-NBARC-LRR or TIR-NBARC-LRR) from 02_nlr_annotation/<genome>/*_nlr.txt
     -> (seqid, start, end) per locus.
  2. Reads ANNEVO's gene coordinates from 00_genomes/<genome>.annevo.gff3
     -> gene_id -> (seqid, start, end).
  3. Reads HMMER's domtblout to get which protein IDs are "complete-like"
     (NB-ARC + at least one LRR domain), and maps protein_id -> gene_id
     by stripping the transcript suffix (e.g. "...-g12.t1" -> "...-g12").
  4. For each genome, checks genomic overlap (same seqid, overlapping
     start-end range) between NLR-Annotator loci and ANNEVO complete-like
     gene loci, and reports how many loci match between the two methods.

Usage:
    cd ~/project
    python3 compare_nlr_loci.py
"""

import re
from pathlib import Path
from collections import defaultdict
import pandas as pd

NLR_ANNOTATOR_DIR = Path("02_nlr_annotation")
ANNEVO_GFF_DIR = Path("00_genomes")
DOMTBL_DIR = Path("nlr_domain_search/per_species")
OUT_CSV = Path("nlr_locus_comparison.csv")

LRR_PROFILES = {"LRR_1", "LRR_3", "LRR_8", "LRR_4"}
COMPLETE_ARCHITECTURES = {"CC-NBARC-LRR", "TIR-NBARC-LRR"}


def load_nlr_annotator_loci(genome):
    """Return list of (seqid, start, end) for complete NLR-Annotator candidates."""
    species_dir = NLR_ANNOTATOR_DIR / genome
    if not species_dir.is_dir():
        return None
    nlr_files = list(species_dir.glob("*_nlr.txt"))
    if not nlr_files:
        return None
    loci = []
    with open(nlr_files[0]) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 5:
                continue
            seqid, _id, arch, start, end = fields[0], fields[1], fields[2], fields[3], fields[4]
            if arch in COMPLETE_ARCHITECTURES:
                loci.append((seqid, int(start), int(end)))
    return loci


def load_annevo_gene_coords(genome):
    """Return dict: gene_id -> (seqid, start, end) from the ANNEVO gff3."""
    gff_path = ANNEVO_GFF_DIR / f"{genome}.annevo.gff3"
    if not gff_path.exists():
        return None
    genes = {}
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            seqid, start, end, attrs = fields[0], fields[3], fields[4], fields[8]
            m = re.search(r"ID=([^;]+)", attrs)
            if not m:
                continue
            gene_id = m.group(1)
            genes[gene_id] = (seqid, int(start), int(end))
    return genes


def load_complete_like_gene_ids(genome):
    """Return set of gene IDs with NB-ARC + LRR in the same protein."""
    domtbl_path = DOMTBL_DIR / f"{genome}.domtblout"
    if not domtbl_path.exists():
        return None
    protein_domains = defaultdict(set)
    with open(domtbl_path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.split()
            hmm_name, protein_id = fields[0], fields[3]
            protein_domains[protein_id].add(hmm_name)

    gene_ids = set()
    for protein_id, domains in protein_domains.items():
        if "NB-ARC" in domains and (domains & LRR_PROFILES):
            gene_id = re.sub(r"\.t\d+$", "", protein_id)
            gene_ids.add(gene_id)
    return gene_ids


def overlaps(a_start, a_end, b_start, b_end):
    return a_start <= b_end and b_start <= a_end


def match_loci(nlr_loci, annevo_loci):
    """Return (n_nlr_matched, n_annevo_matched) via genomic overlap."""
    by_seqid = defaultdict(list)
    for seqid, start, end in annevo_loci:
        by_seqid[seqid].append((start, end))
    for seqid in by_seqid:
        by_seqid[seqid].sort()

    nlr_matched = 0
    annevo_matched_set = set()

    for seqid, start, end in nlr_loci:
        candidates = by_seqid.get(seqid, [])
        found = False
        for i, (a_start, a_end) in enumerate(candidates):
            if overlaps(start, end, a_start, a_end):
                found = True
                annevo_matched_set.add((seqid, a_start, a_end))
        if found:
            nlr_matched += 1

    return nlr_matched, len(annevo_matched_set)


rows = []
genomes = sorted(p.name for p in NLR_ANNOTATOR_DIR.iterdir() if p.is_dir())

print(f"Found {len(genomes)} species folders in {NLR_ANNOTATOR_DIR}.\n")

for genome in genomes:
    nlr_loci = load_nlr_annotator_loci(genome)
    gene_coords = load_annevo_gene_coords(genome)
    complete_like_ids = load_complete_like_gene_ids(genome)

    if nlr_loci is None or gene_coords is None or complete_like_ids is None:
        print(f"  [!] Skipping {genome}: missing input file(s).")
        continue

    annevo_loci = [gene_coords[gid] for gid in complete_like_ids if gid in gene_coords]
    unmapped = len(complete_like_ids) - len(annevo_loci)

    n_nlr_matched, n_annevo_matched = match_loci(nlr_loci, annevo_loci)

    n_nlr_total = len(nlr_loci)
    n_annevo_total = len(annevo_loci)

    rows.append({
        "genome": genome,
        "nlr_annotator_complete_loci": n_nlr_total,
        "annevo_complete_like_loci": n_annevo_total,
        "annevo_unmapped_gene_ids": unmapped,
        "nlr_annotator_loci_matched": n_nlr_matched,
        "annevo_loci_matched": n_annevo_matched,
        "nlr_annotator_loci_unmatched": n_nlr_total - n_nlr_matched,
        "pct_nlr_annotator_confirmed_by_annevo": round(100 * n_nlr_matched / n_nlr_total, 1) if n_nlr_total else None,
        "pct_annevo_confirmed_by_nlr_annotator": round(100 * n_annevo_matched / n_annevo_total, 1) if n_annevo_total else None,
    })
    print(f"  {genome}: {n_nlr_matched}/{n_nlr_total} NLR-Annotator loci confirmed by ANNEVO+HMMER "
          f"({rows[-1]['pct_nlr_annotator_confirmed_by_annevo']}%)")

df = pd.DataFrame(rows).sort_values("genome")
df.to_csv(OUT_CSV, index=False)

print()
print("=" * 70)
total_nlr = df["nlr_annotator_complete_loci"].sum()
total_matched = df["nlr_annotator_loci_matched"].sum()
total_annevo = df["annevo_complete_like_loci"].sum()
total_annevo_matched = df["annevo_loci_matched"].sum()

print(f"TOTAL across {len(df)} species:")
print(f"  NLR-Annotator complete loci:                 {total_nlr}")
print(f"  ...confirmed by ANNEVO+HMMER (overlap):       {total_matched} ({100*total_matched/total_nlr:.1f}%)")
print(f"  ANNEVO+HMMER complete-like loci:              {total_annevo}")
print(f"  ...confirmed by NLR-Annotator (overlap):      {total_annevo_matched} ({100*total_annevo_matched/total_annevo:.1f}%)")
print()
print(f"Saved per-species results to: {OUT_CSV}")
print("=" * 70)
