#!/usr/bin/env python3
"""Summarise the HyPhy FEL/MEME runs for the selected orthogroups (Results 4.3).

For each orthogroup it reports: copies in the original tree and analysed after
codon recovery, genomes, focal species (largest contributor in the original
tree), its ploidy and share of the orthogroup, polyploid share of the analysed
copies, alignment size (codons, empty codons, informative codons with <50% gaps)
and the number of sites called by FEL (positive/negative) and MEME at p <= 0.1,
both raw and restricted to informative codons.

Run from ~/project:
    python3 selection_summary.py --out selection_summary
Outputs <out>.csv (one row per orthogroup) and <out>_composition.csv
(one row per orthogroup x species).
"""
import argparse, csv, json, re
from collections import Counter, defaultdict

OGS = ["OG0000000", "OG0000029", "OG0000018", "OG0000008", "OG0000038"]
P_THRESHOLD = 0.1          # HyPhy default for FEL and MEME
GAP_MAX = 0.5              # informative codon: < 50% gapped sequences
NAME = re.compile(r"^(GC[AF]_\d+_\d)_(.+?)_genomic_[^|]+\|(.+?)_nlr\d+\|")


def parse(name):
    m = NAME.match(name)
    if not m:
        return name[:25], "?"
    return m.group(2).replace("_", " "), m.group(3)


def col(headers, key):
    for i, h in enumerate(headers):
        if h[0] == key:
            return i
    raise KeyError(key)


def read_fasta(path):
    names, seqs = [], []
    for line in open(path):
        if line.startswith(">"):
            names.append(line[1:].strip()); seqs.append([])
        else:
            seqs[-1].append(line.strip())
    return names, ["".join(s).upper() for s in seqs]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hyphy-dir", default="05_selection_hyphy")
    ap.add_argument("--phylo-dir", default="04_phylogeny")
    ap.add_argument("--ploidy", default="ploidy_binary.csv")
    ap.add_argument("--out", default="selection_summary")
    a = ap.parse_args()

    ploidy = {r[0].strip(): r[1].strip() for r in csv.reader(open(a.ploidy)) if len(r) >= 2}
    rows, comp = [], []
    for og in OGS:
        d = f"{a.hyphy_dir}/{og}"
        names, seqs = read_fasta(f"{d}/{og}_codon_alignment.fasta")
        tips = re.findall(r"[(,]([^(),:;]+):", open(f"{a.phylo_dir}/{og}/{og}.treefile").read())
        orig = Counter(parse(t)[0] for t in tips)
        kept, scaf = Counter(), defaultdict(Counter)
        for n in names:
            s, c = parse(n); kept[s] += 1; scaf[s][c] += 1
        focal, focal_n = orig.most_common(1)[0]
        poly = sum(k for s, k in kept.items() if ploidy.get(s) == "polyploid")
        for s in orig:
            top = scaf[s].most_common(1)[0][1] if kept[s] else 0
            comp.append([og, s, ploidy.get(s, "?"), orig[s], kept[s], len(scaf[s]),
                         round(100 * top / kept[s], 1) if kept[s] else ""])

        n = len(seqs)
        gaps = [sum(1 for s in seqs if s[c:c + 3] == "---") / n for c in range(0, len(seqs[0]), 3)]
        inf = {i + 1 for i, g in enumerate(gaps) if g < GAP_MAX}

        fel = json.load(open(f"{d}/{og}_FEL.json")); meme = json.load(open(f"{d}/{og}_MEME.json"))
        h = fel["MLE"]["headers"]; r = fel["MLE"]["content"]["0"]
        ia, ib, ip = col(h, "alpha"), col(h, "beta"), col(h, "p-value")
        pos = {i + 1 for i, x in enumerate(r) if x[ip] <= P_THRESHOLD and x[ib] > x[ia]}
        neg = {i + 1 for i, x in enumerate(r) if x[ip] <= P_THRESHOLD and x[ib] < x[ia]}
        hm = meme["MLE"]["headers"]; rm = meme["MLE"]["content"]["0"]; ipm = col(hm, "p-value")
        ep = {i + 1 for i, x in enumerate(rm) if x[ipm] <= P_THRESHOLD}
        k = len(inf)
        rows.append({
            "orthogroup": og, "copies_original": len(tips), "copies_analysed": n,
            "genomes_analysed": len(kept), "focal_species": focal,
            "focal_ploidy": ploidy.get(focal, "?"), "focal_copies_original": focal_n,
            "focal_share_pct": round(100 * focal_n / len(tips), 1),
            "polyploid_copies_pct": round(100 * poly / n, 1),
            "codons": len(gaps), "empty_codons": sum(g == 1 for g in gaps), "informative_codons": k,
            "FEL_positive": len(pos), "FEL_negative": len(neg), "MEME": len(ep),
            "shared_FEL_MEME": len(pos & ep),
            "FEL_positive_informative": len(pos & inf), "MEME_informative": len(ep & inf),
            "FEL_positive_pct_informative": round(100 * len(pos & inf) / k, 1),
            "MEME_pct_informative": round(100 * len(ep & inf) / k, 1),
        })

    with open(a.out + ".csv", "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
    with open(a.out + "_composition.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["orthogroup", "species", "ploidy", "copies_original", "copies_analysed",
                    "scaffolds", "top_scaffold_pct"])
        w.writerows(comp)
    for r in rows:
        print(r)


if __name__ == "__main__":
    main()
