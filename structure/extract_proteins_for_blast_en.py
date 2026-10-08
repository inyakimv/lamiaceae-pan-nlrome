#!/usr/bin/env python3
"""
Extracts UNGAPPED protein sequences (removes alignment gaps)
for one or more orthogroups, ready to use as BLAST queries
against RefPlantNLR.

Usage:
    cd ~/project
    python3 extract_proteins_for_blast_en.py OG0000038 OG0000018 OG0000029 OG0000008

Reads from: 03_orthofinder/proteomes/OrthoFinder/Results_*/MultipleSequenceAlignments/<OG>.fa
Writes to:  10_validation/<OG>_proteins_for_blast.fasta
"""
import sys
import os
import glob
from Bio import SeqIO

PROJECT_DIR = os.path.expanduser("~/project")
OUT_DIR = os.path.join(PROJECT_DIR, "10_validation")


def find_msa_file(og_id):
    matches = glob.glob(os.path.join(
        PROJECT_DIR, "03_orthofinder/proteomes/OrthoFinder/Results_*",
        "MultipleSequenceAlignments", f"{og_id}.fa"
    ))
    return matches[0] if matches else None


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 extract_proteins_for_blast_en.py OG0000038 [OG0000018 ...]")
        sys.exit(1)

    os.makedirs(OUT_DIR, exist_ok=True)

    for og_id in sys.argv[1:]:
        msa_file = find_msa_file(og_id)
        if not msa_file:
            print(f"[!] Cannot find alignment for {og_id}, skipping.")
            continue

        out_path = os.path.join(OUT_DIR, f"{og_id}_proteins_for_blast.fasta")
        n = 0
        with open(out_path, "w") as out:
            for rec in SeqIO.parse(msa_file, "fasta"):
                ungapped = str(rec.seq).replace("-", "")
                out.write(f">{og_id}|{rec.id}\n{ungapped}\n")
                n += 1
        print(f"{og_id}: {n} sequences written to {out_path}")


if __name__ == "__main__":
    main()
