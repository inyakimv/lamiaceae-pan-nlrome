#!/usr/bin/env python3
"""Partition the NB-ARC tree (Lamiaceae + Arabidopsis + K. nitens + RefPlantNLR anchors)
into maximal subtrees of <= CAP leaves and annotate each clade with its reference NLRs.

Usage: partition_refs.py TREE TableS1.xlsx OUTPREFIX [CAP]
"""
import sys, re, csv
from collections import Counter
from Bio import Phylo
import openpyxl

tree_path, s1_path, prefix = sys.argv[1:4]
CAP = int(sys.argv[4]) if len(sys.argv) > 4 else 300
sys.setrecursionlimit(20000)


def sanitize(s):
    return re.sub(r"[^A-Za-z0-9_.\-]", "_", s)


# ---- Table S1 -> metadata keyed by the sanitized ID used in the tree
meta = {}
ws = openpyxl.load_workbook(s1_path, read_only=True)["RefNLRs"]
rows = ws.iter_rows(values_only=True)
hdr = next(rows)
ix = {h: i for i, h in enumerate(hdr) if h}
for r in rows:
    if not r[0]:
        continue
    sub = str(r[ix["Subclass"]])
    tok = [x.strip() for x in sub.split(";")]
    # Table S1 group codes: GT = TIR-NLR, GR = RPW8-NLR, G10 / *CNL* = CC-NLR
    cls = ("TNL" if "GT" in tok else "RNL" if "GR" in tok
           else "CNL" if (any("CNL" in x for x in tok) or "G10" in tok) else "other")
    meta[sanitize(str(r[0]))] = dict(id=str(r[0]), subclass=sub, cls=cls, host=str(r[ix["Host"]]),
                                     pathogen=str(r[ix["Pathogen"]]), pathogen_spp=str(r[ix["Pathogen_spp"]]))


def kind(name):
    if name.startswith("REF|"): return "REF"
    sp = name.split("|")[0].split("__")[0]
    if sp.startswith("Athaliana"): return "ATH"
    if sp.startswith("Knitens"): return "ALGA"
    return "LAM"


def genus(name):
    sp = name.split("|")[0].split("__")[0]
    m = re.match(r"GC[AF]_\d+\.\d+_([A-Za-z]+)_", sp)
    return m.group(1) if m else "NA"


t = Phylo.read(tree_path, "newick")
t.root_at_midpoint()
n = {}
for c in t.find_clades(order="postorder"):
    n[c] = 1 if c.is_terminal() else sum(n[k] for k in c.clades)
print(f"{tree_path}: {n[t.root]} leaves; REF with metadata: "
      f"{sum(1 for l in t.get_terminals() if l.name.startswith('REF|') and l.name[4:] in meta)} of "
      f"{sum(1 for l in t.get_terminals() if l.name.startswith('REF|'))}")


def partition(cap):
    units, stack = [], [t.root]
    while stack:
        c = stack.pop()
        if n[c] <= cap: units.append(c)
        else: stack.extend(c.clades)
    return units


def pat(units, k):
    return sorted((sum(1 for l in u.get_terminals() if kind(l.name) == k) for u in units), reverse=True)


for cap in (200, 300, 400):
    u = partition(cap)
    print(f"  cap {cap}: {len(u)} units; Arabidopsis per unit: {[x for x in pat(u, 'ATH') if x]}")

units = partition(CAP)
units.sort(key=lambda c: -n[c])
info, members, anchors, collapsed = [], [], [], {}
for i, u in enumerate(units, 1):
    cid = f"Clade{i:02d}"
    leaves = u.get_terminals()
    kinds = Counter(kind(l.name) for l in leaves)
    g = Counter(genus(l.name) for l in leaves if kind(l.name) == "LAM")
    dom, dn = (g.most_common(1)[0] if g else ("NA", 0))
    refs = [l.name[4:] for l in leaves if kind(l.name) == "REF"]
    rc = Counter(meta[r]["cls"] for r in refs if r in meta)
    info.append(dict(clade_id=cid, size=len(leaves), n_lamiaceae=kinds["LAM"], dominant_genus=dom,
                     dominant_pct=round(100 * dn / max(1, kinds["LAM"])),
                     n_arabidopsis=kinds["ATH"], n_alga=kinds["ALGA"], n_ref=kinds["REF"],
                     ref_classes=";".join(f"{k}:{v}" for k, v in sorted(rc.items())),
                     ref_ids=";".join(sorted(meta[r]["id"] if r in meta else r for r in refs)),
                     stem_support=("" if u.confidence is None else u.confidence)))
    for l in leaves:
        members.append((cid, l.name))
    for r in refs:
        m = meta.get(r, {})
        anchors.append((cid, m.get("id", r), m.get("cls", "?"), m.get("subclass", "?"), m.get("host", "?"),
                        m.get("pathogen", "?"), m.get("pathogen_spp", "?")))
    collapsed[u] = cid


def collapse(c):
    if c in collapsed:
        c.clades = []; c.name = collapsed[c]; c.confidence = None; return
    for k in list(c.clades): collapse(k)


collapse(t.root)
Phylo.write(t, prefix + "_collapsed.treefile", "newick")
with open(prefix + "_clade_info.csv", "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(info[0].keys())); w.writeheader(); w.writerows(info)
with open(prefix + "_clade_members.tsv", "w") as f:
    for cid, name in members: f.write(f"{cid}\t{name}\n")
with open(prefix + "_anchors.tsv", "w") as f:
    f.write("clade\tref_id\tclass\tsubclass\thost\tpathogen_group\tpathogen_species\n")
    for a in sorted(anchors): f.write("\t".join(a) + "\n")
sizes = [r["size"] for r in info]
print(f"CAP {CAP}: {len(info)} clades; max {max(sizes)}, median {sorted(sizes)[len(sizes)//2]}; sum {sum(sizes)}; "
      f"Arabidopsis {sum(r['n_arabidopsis'] for r in info)}; alga {sum(r['n_alga'] for r in info)}; "
      f"refs {sum(r['n_ref'] for r in info)}; Lamiaceae {sum(r['n_lamiaceae'] for r in info)}")
print("Refs by class:", Counter(a[2] for a in anchors))
