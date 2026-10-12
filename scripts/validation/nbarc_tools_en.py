#!/usr/bin/env python3
"""
Helper for the NB-ARC phylogeny of the 478 Lamiaceae proteins plus RefPlantNLR
(Sections 3.8 / 4.4, Annex A20). Four sub-commands, called by
run_nbarc_phylogeny_en.sh:

  prepare  join query proteins and RefPlantNLR, give short safe IDs, write idmap.tsv
  extract  keep the NB-ARC domain (best hmmscan hit) of each protein
  trim     remove alignment columns that are mostly gaps (and very incomplete sequences)
  placement  position of every query in the tree: smallest clade with RefPlantNLR
           proteins, its references and bootstrap support (topology only)

Requires: Python 3, Biopython (placement only).
"""
import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict


# ---------------------------------------------------------------- helpers
def read_fasta(path):
    name, seq = None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        yield name, "".join(seq)


def write_fasta(path, records, width=80):
    with open(path, "w") as out:
        for name, seq in records:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def reference_group(name):
    """Same rule as summarize_blast_results_en.py."""
    if name.startswith("Rpi-amr1"):
        return "Rpi-amr1"
    if name.startswith("Rpi-chc"):
        return "Rpi-chc1"
    if name.startswith("NRC0"):
        return "NRC0"
    if "NRC" in name:
        return "other NRC"
    return "other"


# ---------------------------------------------------------------- prepare
def cmd_prepare(a):
    os.makedirs(a.outdir, exist_ok=True)
    recs, idmap = [], []
    n_q = n_r = 0
    for path in sorted(glob.glob(a.queries)):
        for header, seq in read_fasta(path):
            n_q += 1
            new = f"Q{n_q:04d}"
            orig = header.split()[0]               # same as BLAST qseqid
            og = orig.split("|")[0]
            recs.append((new, seq.replace("*", "").replace("-", "")))
            idmap.append((new, "query", og, orig, og))
    for header, seq in read_fasta(a.ref):
        n_r += 1
        new = f"R{n_r:04d}"
        label = header.split()[0]
        recs.append((new, seq.replace("*", "").replace("-", "")))
        idmap.append((new, "reference", reference_group(header), label, header))
    write_fasta(os.path.join(a.outdir, "all_proteins.fasta"), recs)
    with open(os.path.join(a.outdir, "idmap.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "source", "group_or_orthogroup", "original_id", "original_header"])
        w.writerows(idmap)
    print(f"prepare: {n_q} query proteins + {n_r} reference proteins -> {a.outdir}/all_proteins.fasta")


# ---------------------------------------------------------------- extract
def cmd_extract(a):
    # hmmscan --domtblout: target(profile) first, query(protein) fourth
    best = {}
    with open(a.domtblout) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            c = line.split()
            profile, prot = c[0], c[3]
            if profile != a.profile:
                continue
            score = float(c[13])                    # domain score
            env_from, env_to = int(c[19]), int(c[20])
            if prot not in best or score > best[prot][0]:
                best[prot] = (score, env_from, env_to)
    idmap = {}
    with open(os.path.join(a.workdir, "idmap.tsv")) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            idmap[row["id"]] = row
    kept, rows = [], []
    for name, seq in read_fasta(os.path.join(a.workdir, "all_proteins.fasta")):
        if name in best:
            _, s, e = best[name]
            dom = seq[s - 1:e]
            if len(dom) >= a.min_len:
                kept.append((name, dom))
                rows.append((name, idmap[name]["source"], s, e, len(dom)))
                continue
        rows.append((name, idmap[name]["source"], "", "", 0))
    write_fasta(os.path.join(a.workdir, "nbarc_domains.fasta"), kept)
    with open(os.path.join(a.workdir, "nbarc_domain_coordinates.tsv"), "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["id", "source", "env_from", "env_to", "domain_length"])
        w.writerows(rows)
    c_all = Counter(idmap[n]["source"] for n, _ in read_fasta(os.path.join(a.workdir, "all_proteins.fasta")))
    c_kept = Counter(idmap[n]["source"] for n, _ in kept)
    for src in ("query", "reference"):
        print(f"extract: {src}: {c_kept[src]} of {c_all[src]} proteins with an NB-ARC domain >= {a.min_len} aa")
    lens = sorted(len(s) for _, s in kept)
    if lens:
        print(f"extract: domain length min/median/max = {lens[0]}/{lens[len(lens)//2]}/{lens[-1]}")


# ---------------------------------------------------------------- trim
def cmd_trim(a):
    recs = list(read_fasta(a.aln))
    if not recs:
        sys.exit("ERROR: empty alignment")
    L = len(recs[0][1])
    n = len(recs)
    keep_cols = []
    for j in range(L):
        gaps = sum(1 for _, s in recs if s[j] == "-")
        if gaps / n <= a.max_gap_fraction:
            keep_cols.append(j)
    trimmed = [(nm, "".join(s[j] for j in keep_cols)) for nm, s in recs]
    out, dropped = [], []
    for nm, s in trimmed:
        cov = sum(1 for ch in s if ch != "-") / max(1, len(s))
        (out if cov >= a.min_seq_coverage else dropped).append((nm, s))
    write_fasta(a.out, out)
    print(f"trim: {L} columns -> {len(keep_cols)} (columns with <= {a.max_gap_fraction:.0%} gaps kept)")
    print(f"trim: {len(out)} sequences kept, {len(dropped)} removed (< {a.min_seq_coverage:.0%} of the kept columns)")
    if dropped:
        with open(a.out + ".removed.txt", "w") as f:
            for nm, _ in dropped:
                f.write(nm + "\n")


# ---------------------------------------------------------------- placement
def cmd_placement(a):
    """Position of the query proteins in the NB-ARC tree relative to RefPlantNLR.

    The tree is rooted on the clade that contains the TIR-NLR references (the
    first block of RefPlantNLR, up to --outgroup-last). For every query, the
    smallest clade that contains at least one RefPlantNLR protein is recorded
    (its reference proteins and its ultrafast-bootstrap support), and the
    successive larger clades are listed with the references they add.
    Branch lengths are not used.
    """
    from Bio import Phylo
    tree = Phylo.read(a.tree, "newick")
    idmap, order = {}, []
    with open(a.idmap) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            idmap[row["id"]] = row
            if row["source"] == "reference":
                order.append(row["id"])
    label = {k: v["original_id"] for k, v in idmap.items()}
    ref_label = {v["id"]: v["original_id"] for v in idmap.values() if v["source"] == "reference"}
    last = [i for i in order if ref_label[i] == a.outgroup_last]
    if not last:
        sys.exit(f"ERROR: outgroup reference {a.outgroup_last} not found in idmap")
    out_ids = set(order[:order.index(last[0]) + 1])
    in_tree = {t.name for t in tree.get_terminals()}
    out_leaves = [t for t in tree.get_terminals() if t.name in out_ids]
    outgroup = tree.common_ancestor(out_leaves)
    n_q_out = sum(1 for t in outgroup.get_terminals() if idmap[t.name]["source"] == "query")
    n_out = len(outgroup.get_terminals())
    print(f"placement: outgroup = smallest clade containing {len(out_leaves)} TIR-NLR references "
          f"({n_out} leaves, {n_q_out} query proteins inside)")
    if n_q_out:
        print("WARNING: the outgroup contains query proteins; check the rooting")
    tree.root_with_outgroup(outgroup)

    parent = {}
    for c in tree.find_clades(order="level"):
        for ch in c.clades:
            parent[id(ch)] = c

    def leaf_names(c):
        return [t.name for t in c.get_terminals()]

    blast = {}
    with open(a.blast_summary) as f:
        for row in csv.DictReader(f):
            blast[row["query"]] = row["reference_group"]

    clades = {}                          # id(node) -> info
    per_query = []
    for leaf in tree.get_terminals():
        if idmap[leaf.name]["source"] != "query":
            continue
        node = leaf
        while True:
            node = parent[id(node)]
            names = leaf_names(node)
            refs = [n for n in names if idmap[n]["source"] == "reference"]
            if refs:
                break
        info = clades.setdefault(id(node), {"node": node, "queries": [],
                                            "refs": sorted(ref_label[r] for r in refs)})
        info["queries"].append(leaf.name)
        per_query.append((leaf.name, id(node)))

    ordered = sorted(clades.items(), key=lambda kv: -len(kv[1]["queries"]))
    cid = {k: f"C{i + 1}" for i, (k, _) in enumerate(ordered)}
    rows = []
    for q, k in sorted(per_query):
        info = clades[k]
        refs = info["refs"]
        rows.append({
            "orthogroup": idmap[q]["group_or_orthogroup"], "query": idmap[q]["original_id"],
            "clade": cid[k],
            "clade_support": info["node"].confidence if info["node"].confidence is not None else "",
            "clade_queries": len(info["queries"]), "clade_references": len(refs),
            "clade_reference_names": ";".join(refs[:12]) + (";..." if len(refs) > 12 else ""),
            "blast_group": blast.get(idmap[q]["original_id"], ""),
        })
    with open(a.out, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f"placement: {len(rows)} query proteins in the tree -> {a.out}")

    print("\nSmallest clade containing reference proteins, by clade (C1 = largest):")
    for k, info in ordered:
        node = info["node"]
        og = Counter(idmap[q]["group_or_orthogroup"] for q in info["queries"])
        bl = Counter(blast.get(idmap[q]["original_id"], "?") for q in info["queries"])
        refs = info["refs"]
        print(f"\n  {cid[k]}: support {node.confidence}; {len(info['queries'])} query proteins "
              f"{dict(sorted(og.items()))}")
        print(f"      references ({len(refs)}): " + ", ".join(refs[:14]) + (", ..." if len(refs) > 14 else ""))
        print(f"      BLAST best-hit groups of these queries: {dict(bl.most_common())}")
        seen = set(r for r in leaf_names(node) if idmap[r]["source"] == "reference")
        anc, steps = node, 0
        while id(anc) in parent and steps < a.steps:
            anc = parent[id(anc)]
            now = set(r for r in leaf_names(anc) if idmap[r]["source"] == "reference")
            new = sorted(ref_label[r] for r in now - seen)
            if new:
                steps += 1
                print(f"      + clade [{anc.confidence}] adds {len(new)} references: "
                      + ", ".join(new[:10]) + (", ..." if len(new) > 10 else ""))
                seen = now


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sp = p.add_subparsers(dest="cmd", required=True)
    s = sp.add_parser("prepare")
    s.add_argument("--queries", required=True, help="glob, e.g. '10_validation/OG*_proteins_for_blast.fasta'")
    s.add_argument("--ref", required=True)
    s.add_argument("--outdir", required=True)
    s.set_defaults(f=cmd_prepare)
    s = sp.add_parser("extract")
    s.add_argument("--workdir", required=True)
    s.add_argument("--domtblout", required=True)
    s.add_argument("--profile", default="NB-ARC")
    s.add_argument("--min-len", type=int, default=150)
    s.set_defaults(f=cmd_extract)
    s = sp.add_parser("trim")
    s.add_argument("--aln", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--max-gap-fraction", type=float, default=0.5)
    s.add_argument("--min-seq-coverage", type=float, default=0.5)
    s.set_defaults(f=cmd_trim)
    s = sp.add_parser("placement")
    s.add_argument("--tree", required=True)
    s.add_argument("--idmap", required=True)
    s.add_argument("--blast-summary", required=True)
    s.add_argument("--out", required=True)
    s.add_argument("--outgroup-last", default="TNL2-PN40024",
                   help="last RefPlantNLR entry of the TIR-NLR block used as outgroup")
    s.add_argument("--steps", type=int, default=4, help="number of larger clades listed per clade")
    s.set_defaults(f=cmd_placement)
    a = p.parse_args()
    a.f(a)


if __name__ == "__main__":
    main()
