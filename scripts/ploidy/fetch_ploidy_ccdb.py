#!/usr/bin/env python3
"""
Look up published chromosome counts for each species in the Chromosome Counts
Database (CCDB, Rice et al. 2015; http://ccdb.tau.ac.il) and propose a ploidy level.

It does NOT decide ploidy for you. Species with conflicting counts, no counts, or a
non-integer 2n / x are flagged for manual review. Check them against the genome paper
of each assembly (the ploidy of the sequenced individual is what counts).

Input
    --species        CSV with a column "species"
    --base-numbers   optional CSV with columns: species, x   and/or   genus, x
    --aliases        optional CSV with columns: species, query_name
Output (prefix --out)
    <out>_summary.csv, <out>_raw/<species>.csv, <out>_x_template.csv
Zero or missing counts are ignored.
"""
import argparse
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

import pandas as pd

BASE = "http://ccdb.tau.ac.il/services/countsFull/"


def query_names(species, alias=None):
    names = ([alias] if alias else []) + [species]
    if re.search(r"\sx\s", species):
        names.append(re.sub(r"\sx\s", " ", species))
        names.append(re.sub(r"\sx\s", " × ", species))
    seen, out = set(), []
    for n in names:
        if n not in seen:
            seen.add(n)
            out.append(n)
    return out


def download(species, raw_dir, alias=None, pause=1.0):
    path = raw_dir / (species.replace(" ", "_") + ".csv")
    if path.exists() and path.stat().st_size > 0:
        return path, "cached"
    for name in query_names(species, alias):
        url = BASE + "?" + urllib.parse.urlencode({"species": name.replace(" ", "_"), "format": "csv"})
        try:
            with urllib.request.urlopen(url, timeout=60) as r:
                body = r.read().decode("utf-8", errors="replace")
        except Exception as e:
            print(f"  [!] {species}: request failed ({e})", file=sys.stderr)
            continue
        finally:
            time.sleep(pause)
        if body.strip() and len(body.strip().splitlines()) > 1:
            path.write_text(body)
            return path, f"downloaded as '{name}'"
    path.write_text("")
    return path, "no records"


def find_col(df, pattern, explicit=None):
    if explicit:
        return explicit if explicit in df.columns else None
    hits = [c for c in df.columns if re.search(pattern, str(c), re.I)]
    return hits[0] if hits else None


def to_numbers(series):
    v = pd.to_numeric(series, errors="coerce").dropna()
    return v[v > 0]


def summarise(df, n_col, s_col):
    vals = []
    if s_col:
        vals += [int(round(v)) for v in to_numbers(df[s_col])]
    if n_col:
        vals += [int(round(2 * v)) for v in to_numbers(df[n_col])]
    return Counter(vals)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--species", required=True)
    ap.add_argument("--base-numbers", default=None)
    ap.add_argument("--aliases", default=None)
    ap.add_argument("--out", default="ccdb")
    ap.add_argument("--n-col", default=None)
    ap.add_argument("--2n-col", dest="s_col_arg", default=None)
    args = ap.parse_args()

    raw_dir = Path(args.out + "_raw")
    raw_dir.mkdir(parents=True, exist_ok=True)

    species = pd.read_csv(args.species)["species"].dropna().astype(str).str.strip().tolist()
    base_genus, base_species = {}, {}
    if args.base_numbers:
        b = pd.read_csv(args.base_numbers)
        for _, r in b.iterrows():
            if pd.isna(r.get("x")):
                continue
            if "species" in b.columns and isinstance(r.get("species"), str):
                base_species[r["species"].strip().lower()] = float(r["x"])
            elif "genus" in b.columns and isinstance(r.get("genus"), str):
                base_genus[r["genus"].strip().lower()] = float(r["x"])
    aliases = {}
    if args.aliases:
        a = pd.read_csv(args.aliases)
        aliases = {str(sp).strip(): str(q).strip() for sp, q in zip(a["species"], a["query_name"])}

    rows = []
    for sp in species:
        path, status = download(sp, raw_dir, aliases.get(sp))
        print(f"{sp}: {status}")
        row = {"species": sp, "n_records": 0, "observed_2n": "", "mode_2n": "",
               "x": "", "ploidy_proposed": "", "flag": ""}
        if status == "no records":
            row["flag"] = "no CCDB records - look up the genome paper"
            rows.append(row)
            continue
        df = pd.read_csv(path)
        n_col = find_col(df, r"gameto", args.n_col)
        s_col = find_col(df, r"sporo", args.s_col_arg)
        if not n_col and not s_col:
            sys.exit("Could not find count columns. Columns received from CCDB:\n  "
                     + ", ".join(map(str, df.columns))
                     + "\nRe-run with --n-col NAME and/or --2n-col NAME.")
        counts = summarise(df, n_col, s_col)
        row["n_records"] = len(df)
        if not counts:
            row["flag"] = "records found but no numeric counts"
            rows.append(row)
            continue
        row["observed_2n"] = "; ".join(f"{k} (x{v})" for k, v in sorted(counts.items()))
        mode, _ = counts.most_common(1)[0]
        row["mode_2n"] = mode
        flags = []
        if len(counts) > 1:
            flags.append("several 2n values - check which applies to the sequenced individual")
        genus = sp.split()[0].lower()
        x = base_species.get(sp.lower(), base_genus.get(genus))
        if x is not None:
            row["x"] = x
            ratio = mode / x
            row["ploidy_proposed"] = f"{int(round(ratio))}x"
            if abs(ratio - round(ratio)) > 0.01 or round(ratio) < 1:
                flags.append(f"2n/x = {ratio:.2f} is not an integer - check x or dysploidy")
        elif args.base_numbers:
            flags.append("no base number x given for this species or genus")
        row["flag"] = " | ".join(flags)
        rows.append(row)

    out = pd.DataFrame(rows)
    out.to_csv(args.out + "_summary.csv", index=False)
    out[["species", "mode_2n"]].assign(x="").to_csv(args.out + "_x_template.csv", index=False)
    print(f"\nWrote {args.out}_summary.csv ({len(out)} species), {args.out}_x_template.csv "
          f"and raw records in {raw_dir}/")
    flagged = out[out["flag"] != ""]
    print(f"{len(flagged)} species need a manual check (see the 'flag' column).")


if __name__ == "__main__":
    main()
