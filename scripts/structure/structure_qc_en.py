#!/usr/bin/env python3
"""
Quality control for AlphaFold-predicted structures: reads the
per-residue pLDDT confidence (from the .cif file's B-factor
column) and, if available, the PAE matrix (from the
confidences.json file), and produces a summary CSV classifying
each structure/domain as high/medium/low confidence.

Expects one subfolder per candidate inside
06_structure_retrieval/, each containing the files downloaded
from AlphaFold Server (a .cif model file and, ideally, a
confidences.json file).

Usage:
    cd ~/project
    python3 structure_qc_en.py

pLDDT interpretation (standard AlphaFold thresholds):
    >= 90        : very high confidence
    70-90        : confident
    50-70        : low confidence
    < 50         : very low confidence (likely disordered)
"""
import os
import glob
import json
import csv
import warnings

from Bio.PDB import MMCIFParser

PROJECT_DIR = os.path.expanduser("~/project")
STRUCTURE_DIR = os.path.join(PROJECT_DIR, "06_structure_retrieval")


def find_best_model_index(folder):
    """AlphaFold Server generates 5 predictions per job (model_0
    to model_4), each with a matching summary_confidences_N.json
    that holds a 'ranking_score'. This picks the index with the
    highest ranking score (falls back to index 0 if scores are
    unavailable)."""
    summary_files = glob.glob(os.path.join(folder, "*summary_confidences_*.json"))
    if not summary_files:
        return 0

    best_index, best_score = 0, float("-inf")
    for path in summary_files:
        base = os.path.basename(path)
        try:
            idx = int(base.rsplit("_", 1)[-1].replace(".json", ""))
        except ValueError:
            continue
        try:
            with open(path) as f:
                data = json.load(f)
            score = data.get("ranking_score", 0)
        except Exception:
            score = 0
        if score > best_score:
            best_score = score
            best_index = idx
    return best_index


def find_cif_file(folder, index):
    matches = glob.glob(os.path.join(folder, f"*model_{index}.cif"))
    if matches:
        return matches[0]
    # fallback: any .cif file
    matches = glob.glob(os.path.join(folder, "*.cif"))
    return matches[0] if matches else None


def find_confidences_json(folder, index):
    """The per-residue PAE data is in the 'full_data_N.json' file
    (not a file simply named 'confidences.json')."""
    matches = glob.glob(os.path.join(folder, f"*full_data_{index}.json"))
    if matches:
        return matches[0]
    matches = glob.glob(os.path.join(folder, "*full_data*.json"))
    return matches[0] if matches else None


def get_plddt_per_residue(cif_path):
    """Extracts per-residue pLDDT from the CA atom's B-factor
    column of the CIF file (standard AlphaFold convention)."""
    parser = MMCIFParser(QUIET=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        structure = parser.get_structure("model", cif_path)

    plddts = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if "CA" in residue:
                    plddts.append(residue["CA"].get_bfactor())
        break  # only the first model
    return plddts


def get_mean_pae(json_path):
    """Tries common key names used across AlphaFold Server/AF3
    confidence JSON files to find the PAE matrix and returns its
    overall mean."""
    if not json_path or not os.path.exists(json_path):
        return None
    with open(json_path) as f:
        data = json.load(f)

    for key in ("pae", "predicted_aligned_error", "pae_matrix"):
        if key in data:
            matrix = data[key]
            flat = [v for row in matrix for v in row]
            return sum(flat) / len(flat) if flat else None
    return None


def classify_plddt(mean_plddt):
    if mean_plddt >= 90:
        return "very_high"
    elif mean_plddt >= 70:
        return "confident"
    elif mean_plddt >= 50:
        return "low"
    else:
        return "very_low"


def main():
    if not os.path.isdir(STRUCTURE_DIR):
        print(f"ERROR: cannot find {STRUCTURE_DIR}")
        return

    subfolders = [f for f in glob.glob(os.path.join(STRUCTURE_DIR, "*")) if os.path.isdir(f)]
    if not subfolders:
        print(f"ERROR: no subfolders found inside {STRUCTURE_DIR}")
        print("Extract each AlphaFold Server .zip into its own subfolder first.")
        return

    results = []
    for folder in sorted(subfolders):
        name = os.path.basename(folder)
        best_index = find_best_model_index(folder)
        cif_path = find_cif_file(folder, best_index)
        if not cif_path:
            print(f"[!] {name}: no .cif file found, skipping.")
            continue

        try:
            plddts = get_plddt_per_residue(cif_path)
        except Exception as e:
            print(f"[!] {name}: error reading structure ({e}), skipping.")
            continue

        if not plddts:
            print(f"[!] {name}: no CA atoms found, skipping.")
            continue

        mean_plddt = sum(plddts) / len(plddts)
        n_low = sum(1 for p in plddts if p < 50)
        pct_low = n_low / len(plddts) * 100

        json_path = find_confidences_json(folder, best_index)
        mean_pae = get_mean_pae(json_path)

        classification = classify_plddt(mean_plddt)

        results.append({
            "candidate": name,
            "n_residues": len(plddts),
            "mean_plddt": round(mean_plddt, 1),
            "pct_residues_low_confidence": round(pct_low, 1),
            "mean_pae": round(mean_pae, 1) if mean_pae is not None else "",
            "classification": classification,
        })
        print(f"{name}: mean pLDDT={mean_plddt:.1f} ({classification}), "
              f"{pct_low:.1f}% residues < 50, mean PAE={mean_pae}")

    out_csv = os.path.join(STRUCTURE_DIR, "structure_qc_summary.csv")
    with open(out_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=[
            "candidate", "n_residues", "mean_plddt",
            "pct_residues_low_confidence", "mean_pae", "classification"
        ])
        writer.writeheader()
        writer.writerows(results)

    print(f"\nSummary saved to: {out_csv}")
    print(f"Structures processed: {len(results)}")


if __name__ == "__main__":
    main()
