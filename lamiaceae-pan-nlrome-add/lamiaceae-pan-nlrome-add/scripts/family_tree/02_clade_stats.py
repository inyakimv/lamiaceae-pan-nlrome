#!/usr/bin/env python3
"""Per-clade conservation statistics for the collapsed NB-ARC tree.

Expands the cd-hit (85%) representatives with their cluster members and reports, per clade,
the number of Lamiaceae species present, genes, and genes per species.

Usage: 02_clade_stats.py CDHIT_CLSTR CLADE_MEMBERS_TSV CLADE_INFO_CSV OUTPREFIX
"""
import re, sys, csv, statistics, collections
import pandas as pd
clstr, members_tsv, info_csv, prefix = sys.argv[1:5]
members, cur = {}, []
def flush():
    global cur
    rep = [m for m, star in cur if star]
    if rep: members[rep[0]] = [m for m, _ in cur]
    cur = []
for line in open(clstr):
    if line.startswith(">Cluster"): flush(); continue
    m = re.search(r">(.+?)\.\.\.", line)
    if m: cur.append((m.group(1), "*" in line.split("...")[1][:6]))
flush()
clade = {}
for line in open(members_tsv):
    c, leaf = line.rstrip("\n").split("\t"); clade[leaf] = c
info = pd.read_csv(info_csv).set_index("clade_id")
species = lambda name: name.split("|")[0].split("_genomic")[0]
per = collections.defaultdict(collections.Counter)
for rep, ms in members.items():
    c = clade.get(rep)
    if not c: continue
    for m in ms:
        if not m.startswith(("Athaliana", "Knitens", "REF")): per[c][species(m)] += 1
rows = []
for c, cnt in per.items():
    cp = list(cnt.values())
    rows.append(dict(clade=c, n_species=len(cnt), n_genes=sum(cp), median_copies=statistics.median(cp),
                     copies_per_species=round(sum(cp) / len(cp), 2), n_arabidopsis=info.loc[c, "n_arabidopsis"],
                     n_ref=info.loc[c, "n_ref"], ref_classes=info.loc[c, "ref_classes"]))
t = pd.DataFrame(rows).sort_values(["n_species", "n_genes"], ascending=[False, True])
t.to_csv(prefix + "_clade_conservation.tsv", sep="\t", index=False)
t[["clade", "n_genes"]].rename(columns={"n_genes": "lamiaceae_genes_expanded"}).to_csv(prefix + "_clade_genes.tsv", sep="\t", index=False)
print(f"{len(t)} clades with Lamiaceae genes; "
      f"{((t.n_species >= 35) & (t.copies_per_species <= 3)).sum()} conserved low-copy clades")
