#!/usr/bin/env python3
"""
Analyzes the most expanded orthogroup for each of the 46
species: NLR type composition (CNL/TNL/other) and evidence of
tandem duplication (genes physically clustered on the same
scaffold, within a threshold distance).

Usage:
    cd ~/project/03_orthofinder/proteomes/OrthoFinder/Results_Jul07/Orthogroups
    python3 analizar_tandem_46especies.py

Requires:
    - Orthogroups.GeneCount.tsv and Orthogroups.tsv (in the current folder)
    - The *_nlr.txt files from NLR-Annotator in ~/project/02_nlr_annotation/

Adjustable parameter:
    TANDEM_THRESHOLD_BP: maximum distance (in bp) between
    consecutive genes of the same group, on the same scaffold,
    to consider them part of the same tandem cluster.
    Default: 500,000 bp (500 kb) -- adjust based on what you
    want to justify in your methodology.
"""

import glob
import os
import csv
import sys

TANDEM_THRESHOLD_BP = 500_000
NLR_DIR = os.path.expanduser("~/project/02_nlr_annotation")
GENECOUNT_TSV = "Orthogroups.GeneCount.tsv"
ORTHOGROUPS_TSV = "Orthogroups.tsv"
OUTPUT_CSV = "tandem_analysis_46species.csv"


def load_genecount(path):
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        species_cols = header[1:-1]  # excludes "Orthogroup" and "Total"
        rows = []
        for line in f:
            rows.append(line.rstrip("\n").split("\t"))
    return species_cols, rows


def find_top_orthogroup_per_species(species_cols, rows):
    """For each species, finds its row (orthogroup) with the
    highest copy number, and the total for that group."""
    results = {}
    for idx, species in enumerate(species_cols, start=1):
        best_og, best_count, best_total = None, -1, None
        for row in rows:
            try:
                count = int(row[idx])
            except (ValueError, IndexError):
                continue
            if count > best_count:
                best_count = count
                best_og = row[0]
                best_total = int(row[-1])
        results[species] = (best_og, best_count, best_total)
    return results


def load_orthogroups_gene_lists(path):
    """Returns dict: (species_idx) -> dict(orthogroup -> gene list)"""
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        data = {i: {} for i in range(1, len(header))}
        for line in f:
            cols = line.rstrip("\n").split("\t")
            og = cols[0]
            for i in range(1, len(cols)):
                if cols[i].strip():
                    genes = [g.strip() for g in cols[i].split(",") if g.strip()]
                    data[i][og] = genes
    return header, data


def extract_locus_id(gene_field):
    """From 'GCA_XXXX|scaffold_nlrN|method' extracts 'scaffold_nlrN'"""
    parts = gene_field.split("|")
    if len(parts) >= 2:
        return parts[1]
    return gene_field


def load_nlr_annotations(species_genome_name):
    """Finds that species' _nlr.txt file and returns
    dict: locus_id -> (scaffold, type, start, end)"""
    pattern = os.path.join(NLR_DIR, species_genome_name, "*_nlr.txt")
    matches = glob.glob(pattern)
    if not matches:
        return None
    annotations = {}
    with open(matches[0]) as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 5:
                continue
            scaffold, locus_id, nlr_type, start, end = cols[0], cols[1], cols[2], cols[3], cols[4]
            try:
                annotations[locus_id] = (scaffold, nlr_type, int(start), int(end))
            except ValueError:
                continue
    return annotations


def compute_tandem_clusters(gene_coords, threshold_bp):
    """gene_coords: list of (scaffold, start, end).
    Returns: (largest_cluster_size, n_genes_in_any_cluster, n_scaffolds_used)"""
    by_scaffold = {}
    for scaffold, start, end in gene_coords:
        by_scaffold.setdefault(scaffold, []).append((start, end))

    largest_cluster = 1 if gene_coords else 0
    n_in_cluster = 0
    for scaffold, positions in by_scaffold.items():
        positions.sort()
        cluster_size = 1
        for i in range(1, len(positions)):
            gap = positions[i][0] - positions[i - 1][1]
            if gap <= threshold_bp:
                cluster_size += 1
            else:
                if cluster_size >= 2:
                    n_in_cluster += cluster_size
                largest_cluster = max(largest_cluster, cluster_size)
                cluster_size = 1
        if cluster_size >= 2:
            n_in_cluster += cluster_size
        largest_cluster = max(largest_cluster, cluster_size)

    return largest_cluster, n_in_cluster, len(by_scaffold)


def main():
    if not os.path.exists(GENECOUNT_TSV) or not os.path.exists(ORTHOGROUPS_TSV):
        print(f"ERROR: cannot find {GENECOUNT_TSV} or {ORTHOGROUPS_TSV} in the current folder.")
        sys.exit(1)

    species_cols, gc_rows = load_genecount(GENECOUNT_TSV)
    top_og_per_species = find_top_orthogroup_per_species(species_cols, gc_rows)

    og_header, og_gene_lists = load_orthogroups_gene_lists(ORTHOGROUPS_TSV)

    results = []
    for idx, species in enumerate(species_cols, start=1):
        top_og, count_in_species, total_in_group = top_og_per_species[species]
        if top_og is None or count_in_species <= 0:
            continue

        gene_fields = og_gene_lists.get(idx, {}).get(top_og, [])
        loci = [extract_locus_id(g) for g in gene_fields]

        annotations = load_nlr_annotations(species)
        if annotations is None:
            print(f"  [!] Could not find _nlr.txt file for {species}. Skipping.")
            continue

        type_counts = {}
        gene_coords = []
        for locus in loci:
            if locus in annotations:
                scaffold, nlr_type, start, end = annotations[locus]
                type_counts[nlr_type] = type_counts.get(nlr_type, 0) + 1
                gene_coords.append((scaffold, start, end))

        largest_cluster, n_in_cluster, n_scaffolds = compute_tandem_clusters(
            gene_coords, TANDEM_THRESHOLD_BP
        )

        n_genes_found = len(gene_coords)
        pct_in_tandem = (n_in_cluster / n_genes_found * 100) if n_genes_found else 0
        pct_cnl = (type_counts.get("CC-NBARC-LRR", 0) / n_genes_found * 100) if n_genes_found else 0
        pct_tnl = (type_counts.get("TIR-NBARC-LRR", 0) / n_genes_found * 100) if n_genes_found else 0

        results.append({
            "species": species,
            "top_orthogroup": top_og,
            "copies_in_species": count_in_species,
            "total_in_group": total_in_group,
            "pct_concentration": round(count_in_species / total_in_group * 100, 1) if total_in_group else 0,
            "n_scaffolds_involved": n_scaffolds,
            "largest_tandem_cluster": largest_cluster,
            "pct_genes_in_tandem_clusters": round(pct_in_tandem, 1),
            "pct_CNL": round(pct_cnl, 1),
            "pct_TNL": round(pct_tnl, 1),
        })

        print(f"  -> {species}: {top_og} ({count_in_species} copies, "
              f"largest cluster={largest_cluster}, {pct_in_tandem:.1f}% in tandem)")

    with open(OUTPUT_CSV, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "species", "top_orthogroup", "copies_in_species", "total_in_group",
            "pct_concentration", "n_scaffolds_involved", "largest_tandem_cluster",
            "pct_genes_in_tandem_clusters", "pct_CNL", "pct_TNL"
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSummary saved to: {os.path.abspath(OUTPUT_CSV)}")
    print(f"Total species processed: {len(results)}")


if __name__ == "__main__":
    main()
