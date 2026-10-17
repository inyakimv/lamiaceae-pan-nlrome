#!/usr/bin/env python3
"""Build the FASTA for a clade-level tree: all Lamiaceae proteins of the clade (full length) + selected references.

Usage: python3 06_clade_fasta.py CLADE_GENES.tsv PROTEIN_DIR REFS_FASTA REF_REGEX OUT.fasta [OUTGROUP_REGEX]
  REF_REGEX       regex on RefPlantNLR names to include as ingroup references, e.g. 'ADR1|NRG1'
  OUTGROUP_REGEX  optional regex for outgroup references, e.g. '^ZAR1$|^NRC0$'
Leaf names: <Species>|<gene> for Lamiaceae, REF_<name> for references (OUT_ prefix for the outgroup).
"""
import sys, re, glob, os
import pandas as pd
from Bio import SeqIO
genes, pdir, refs, rx, out = sys.argv[1:6]
orx = sys.argv[6] if len(sys.argv) > 6 else None
d = pd.read_csv(genes, sep="\t")
want = {}
for g in d.gene:
    gen, k = g.split("|"); want.setdefault(gen, set()).add(k)
recs = []
for gen, keys in want.items():
    f = f"{pdir}/{gen}.annevo_proteins.fasta"
    sp = "_".join(gen.split("_")[2:]).replace("_genomic", "")
    for r in SeqIO.parse(f, "fasta"):
        if r.id in keys:
            r.id = f"{sp}|{r.id}"; r.description = ""; recs.append(r)
n_l = len(recs)
for r in SeqIO.parse(refs, "fasta"):
    nm = r.id
    if re.search(rx, nm):
        r.id = "REF_" + re.sub(r"\W", "_", nm); r.description = ""; recs.append(r)
    elif orx and re.search(orx, nm):
        r.id = "OUT_" + re.sub(r"\W", "_", nm); r.description = ""; recs.append(r)
SeqIO.write(recs, out, "fasta")
print(f"{n_l} Lamiaceae proteins + {len(recs)-n_l} references -> {out}")
print([r.id for r in recs if not r.id.startswith(tuple(["Ag", "Aj", "Ba", "Ca"])) and (r.id.startswith("REF_") or r.id.startswith("OUT_"))])
