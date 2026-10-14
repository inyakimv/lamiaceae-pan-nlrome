#!/usr/bin/env python3
"""Prune an orthogroup ML tree to the sequences present in its codon alignment.

prepare_hyphy_input_en.py drops copies whose in-frame nucleotide sequence cannot
be recovered; HyPhy needs the same taxa in the tree and the alignment, so the
IQ-TREE tree is pruned to the sequences that remain. Pruning removes the tips and
collapses the nodes left with a single child (branch lengths are summed).

Usage:
    cd ~/project
    python3 prune_tree_to_alignment.py OG0000038
    python3 prune_tree_to_alignment.py OG0000038 --compare   # also check against an existing pruned tree

Reads  04_phylogeny/<OG>/<OG>.treefile and 05_selection_hyphy/<OG>/<OG>_codon_alignment.fasta
Writes 05_selection_hyphy/<OG>/<OG>_pruned.treefile  (or <OG>_pruned_check.treefile with --compare,
       so an existing pruned tree is never overwritten)
Requires Biopython.
"""
import argparse, os, random
from Bio import Phylo, SeqIO


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("og")
    ap.add_argument("--compare", action="store_true",
                    help="write <OG>_pruned_check.treefile and compare it with the existing <OG>_pruned.treefile")
    a = ap.parse_args()
    d = f"05_selection_hyphy/{a.og}"
    keep = {r.id for r in SeqIO.parse(f"{d}/{a.og}_codon_alignment.fasta", "fasta")}
    tree = Phylo.read(f"04_phylogeny/{a.og}/{a.og}.treefile", "newick")
    tips = {t.name for t in tree.get_terminals()}
    missing = keep - tips
    if missing:
        raise SystemExit(f"{len(missing)} alignment sequences are not tips of the tree, e.g. {sorted(missing)[:3]}")
    for name in sorted(tips - keep):
        tree.prune(name)
    out = f"{d}/{a.og}_pruned_check.treefile" if a.compare else f"{d}/{a.og}_pruned.treefile"
    Phylo.write(tree, out, "newick")
    print(f"{a.og}: tree tips {len(tips)} -> {len(keep)} (removed {len(tips) - len(keep)}); written {out}")

    if a.compare:
        ref = Phylo.read(f"{d}/{a.og}_pruned.treefile", "newick")
        rtips = {t.name for t in ref.get_terminals()}
        print("same tip set as existing pruned tree:", rtips == keep)
        if rtips == keep:
            names = sorted(keep); rnd = random.Random(1)
            pairs = [rnd.sample(names, 2) for _ in range(min(300, len(names) * 2))]
            diff = max(abs(tree.distance(x, y) - ref.distance(x, y)) for x, y in pairs)
            print(f"max difference in patristic distance over {len(pairs)} random tip pairs: {diff:.2e}")


if __name__ == "__main__":
    main()
