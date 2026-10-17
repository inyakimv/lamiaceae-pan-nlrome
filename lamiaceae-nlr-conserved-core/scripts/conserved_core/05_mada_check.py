#!/usr/bin/env python3
"""N-terminal (MADA-like) check for the proteins of one clade.

Usage: python3 05_mada_check.py CLADE_GENES.tsv PROTEIN_DIR OUT.tsv [REFS_FASTA]
  CLADE_GENES.tsv  output of 04_inspect_clade.py (column 'gene', e.g. GCA_..._genomic|CM...-g12.t1)
  PROTEIN_DIR      folder with <genome>_genomic.annevo_proteins.fasta
  REFS_FASTA       optional RefPlantNLR protein FASTA; N-termini of ZAR1 and NRC references are printed for comparison
Reports the first 30 residues, an N-terminal motif flag (M[AE][DE][AV] start) and
the MADA-like score: matches to the ZAR1-type MADA consensus MADAVVSFVLEKLGNLLTKEY over the first 21 residues,
allowing the start codon to sit up to 3 residues downstream.
"""
import sys, glob, os, re, collections
import pandas as pd
from Bio import SeqIO

genes, pdir, out = sys.argv[1:4]
refs = sys.argv[4] if len(sys.argv) > 4 else None
CONS = "MADAVVSFVLEKLGNLLTKEY"

def score(nt):
    best = (0, 0)
    for off in range(0, 4):
        s = nt[off:off + len(CONS)]
        m = sum(a == b for a, b in zip(s, CONS))
        if m > best[0]: best = (m, off)
    return best

d = pd.read_csv(genes, sep="\t")
byfile = collections.defaultdict(set)
for g in d.gene:
    byfile[g.split("|")[0]].add(g.split("|")[1])

rows = []; missing = 0
for gen, keys in byfile.items():
    f = f"{pdir}/{gen}.annevo_proteins.fasta"
    if not os.path.exists(f):
        cand = glob.glob(f"{pdir}/{gen}*"); f = cand[0] if cand else None
    if not f: print("no protein file for", gen); missing += len(keys); continue
    found = set()
    for rec in SeqIO.parse(f, "fasta"):
        k = next((k for k in keys if rec.id == k or rec.id.endswith(k) or k in rec.id), None)
        if k:
            s = str(rec.seq); sc, off = score(s[:30])
            rows.append(dict(species=gen.split("_", 2)[2].replace("_genomic", ""), gene=rec.id, length=len(s),
                             nterm30=s[:30], starts_M=s.startswith("M"), mada_like_start=bool(re.match(r"M[AE][DE][AV]", s)),
                             mada_matches_of_21=sc, offset=off))
            found.add(k)
    missing += len(keys - found)
r = pd.DataFrame(rows).sort_values(["mada_matches_of_21"], ascending=False)
r.to_csv(out, sep="\t", index=False)
print(f"{len(r)} proteins found, {missing} genes without a protein match")
print("proteins starting with M[AE][DE][AV]:", int(r.mada_like_start.sum()), "of", len(r))
print("matches to the MADA consensus (of 21): median", r.mada_matches_of_21.median(), "max", r.mada_matches_of_21.max())
print(r.head(12)[["species", "nterm30", "mada_matches_of_21", "length"]].to_string(index=False))
if refs:
    print("\nReference N-termini (ZAR1 / NRC):")
    for rec in SeqIO.parse(refs, "fasta"):
        if re.search(r"ZAR1|NRC", rec.id, re.I):
            s = str(rec.seq); print(f"  {rec.id[:40]:40s} {s[:30]}  matches={score(s[:30])[0]}")
