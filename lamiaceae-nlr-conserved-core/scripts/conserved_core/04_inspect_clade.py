#!/usr/bin/env python3
"""Inspect one NB-ARC clade: member genes (cd-hit clusters expanded), copies per species,
genomic position (contig, gene index), neighbouring clades and reference NLRs in the tree.

Usage: python3 04_inspect_clade.py CLADE PARTITION_DIR CLSTR OUT_PREFIX [NBARC_FASTA]
Example: python3 04_inspect_clade.py Clade76 family_tree_out nbarc_domains_cdhit85.fasta.clstr family_tree_out/Clade76 nbarc_domains_combined.fasta
"""
import sys, re, collections
import pandas as pd
from Bio import Phylo, SeqIO

clade, D, clstr, out = sys.argv[1:5]
fasta = sys.argv[5] if len(sys.argv) > 5 else None

members_of = {}; cur = []
def flush():
    global cur
    reps = [m for m, r in cur if r]
    if reps: members_of[reps[0]] = [m for m, _ in cur]
    cur = []
for l in open(clstr):
    if l.startswith(">Cluster"): flush(); continue
    m = re.search(r">(.+?)\.\.\.", l)
    if m: cur.append((m.group(1), "*" in l.split("...")[1][:6]))
flush()

leaf_clade = {}
for l in open(f"{D}/refs_clade_members.tsv"):
    c, leaf = l.rstrip("\n").split("\t"); leaf_clade[leaf] = c
reps = [x for x, c in leaf_clade.items() if c == clade]

rows = []
for r in reps:
    for g in members_of.get(r, [r]):
        if g.startswith(("REF|", "Athaliana", "Knitens")): continue
        sp = g.split("|")[0].replace("_genomic", "")
        m = re.match(r"(.+)-g(\d+)\.t\d+$", g.split("|")[1])
        rows.append(dict(gene=g, species=sp, contig=m.group(1) if m else "", gene_index=int(m.group(2)) if m else -1, representative=r, is_rep=(g == r)))
d = pd.DataFrame(rows).sort_values(["species", "contig", "gene_index"])
d.to_csv(out + "_genes.tsv", sep="\t", index=False)

per = d.groupby("species").size()
print(f"{clade}: {len(d)} genes, {per.size} species, {per.mean():.2f} copies/species (max {per.max()})")
print("copies per species distribution:", dict(per.value_counts().sort_index()))
print("species with >1 copy on the same contig within 30 gene indices (possible tandem):")
n = 0
for (sp, ct), g in d.groupby(["species", "contig"]):
    idx = sorted(g.gene_index)
    if any(b - a <= 30 for a, b in zip(idx, idx[1:])):
        print(" ", sp, ct, idx); n += 1
print("  ->", n, "contig(s)")
if fasta:
    keep = set(d.gene)
    SeqIO.write((s for s in SeqIO.parse(fasta, "fasta") if s.id in keep), out + "_nbarc.fasta", "fasta")
    print("NB-ARC sequences written:", out + "_nbarc.fasta")

t = Phylo.read(f"{D}/refs_collapsed.treefile", "newick")   # leaves are clade ids
info = pd.read_csv(f"{D}/refs_clade_info.csv").set_index("clade_id")
leaves = {x.name: x for x in t.get_terminals()}
if clade in leaves:
    dist = sorted((t.distance(leaves[clade], x), n) for n, x in leaves.items() if n != clade)[:8]
    print("closest clades in the tree (distance, clade, Lamiaceae leaves, Arabidopsis, reference classes and ids):")
    for dd, n in dist:
        r = info.loc[n]
        print(f"  {dd:.2f} {n} lam={r.n_lamiaceae} ath={r.n_arabidopsis} ref={r.n_ref} {r.ref_classes} {str(r.ref_ids)[:80]}")
    print("nested groups containing this clade (support, clades):")
    for c in reversed(t.get_path(leaves[clade])[:-1][-4:]):
        print(f"  {c.confidence} {[x.name for x in c.get_terminals()]}")
print("reference NLRs inside this clade:", [x for x in reps if x.startswith("REF|")])
