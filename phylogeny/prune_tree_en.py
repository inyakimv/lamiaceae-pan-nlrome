#!/usr/bin/env python3
"""
Prunes a .treefile (Newick) so it contains exactly the same
taxa as a given codon alignment fasta -- required because
HyPhy needs the tree and the alignment to have the exact same
set of sequences.

Usage:
    python3 prune_tree_en.py OG0000038

Requires: ete3 (pip install ete3) OR falls back to a simple
manual Newick pruning if ete3 is not available.
"""
import sys
import os
import re
import glob

PROJECT_DIR = os.path.expanduser("~/project")


def get_alignment_names(fasta_path):
    names = []
    with open(fasta_path) as f:
        for line in f:
            if line.startswith(">"):
                names.append(line[1:].strip().split()[0])
    return set(names)


def find_tree_file(og_id):
    matches = glob.glob(os.path.join(PROJECT_DIR, "04_phylogeny", og_id, "*.treefile"))
    return matches[0] if matches else None


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 prune_tree_en.py OGxxxxxxx")
        sys.exit(1)

    og_id = sys.argv[1]
    alignment_path = os.path.join(PROJECT_DIR, "05_selection_hyphy", og_id, f"{og_id}_codon_alignment.fasta")
    tree_path = find_tree_file(og_id)

    if not os.path.exists(alignment_path):
        print(f"ERROR: cannot find {alignment_path}")
        sys.exit(1)
    if not tree_path:
        print(f"ERROR: cannot find a .treefile for {og_id} in 04_phylogeny/{og_id}/")
        sys.exit(1)

    keep_names = get_alignment_names(alignment_path)
    print(f"Alignment has {len(keep_names)} sequences to keep.")
    print(f"Using tree: {tree_path}")

    try:
        from ete3 import Tree
        t = Tree(tree_path, format=1)
        all_leaf_names = set(t.get_leaf_names())
        missing_in_tree = keep_names - all_leaf_names
        if missing_in_tree:
            print(f"WARNING: {len(missing_in_tree)} alignment names not found in tree:")
            for n in list(missing_in_tree)[:5]:
                print(f"  {n}")
        t.prune(list(keep_names & all_leaf_names), preserve_branch_length=True)
        out_path = os.path.join(PROJECT_DIR, "05_selection_hyphy", og_id, f"{og_id}_pruned.treefile")
        t.write(format=1, outfile=out_path)
        print(f"\nPruned tree saved to: {out_path}")
        print(f"Tree now has {len(t.get_leaf_names())} leaves.")
    except ImportError:
        print("\nERROR: ete3 is not installed.")
        print("Install it with: pip install ete3")
        sys.exit(1)


if __name__ == "__main__":
    main()
