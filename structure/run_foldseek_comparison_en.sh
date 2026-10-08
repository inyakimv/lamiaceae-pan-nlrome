#!/bin/bash
# ==========================================================
# Compares all predicted structures against each other using
# Foldseek, to detect structural similarity even in cases
# where sequence identity is low (as seen in the BLAST vs
# RefPlantNLR results, ~35% average identity).
#
# IMPORTANT: run from ~/project
#
# Requirements (install once, dedicated environment):
#   conda create -n foldseek_env --solver=classic -c bioconda -c conda-forge foldseek -y
#
# Usage:
#   cd ~/project
#   bash run_foldseek_comparison_en.sh
# ==========================================================

set -uo pipefail

STRUCTURE_DIR="06_structure_retrieval"
OUT_DIR="08_structure_comparison"
DB_DIR="${OUT_DIR}/foldseek_db"

mkdir -p "$OUT_DIR"
mkdir -p "$DB_DIR/structures"

if [ ! -d "$STRUCTURE_DIR" ]; then
    echo "ERROR: cannot find $STRUCTURE_DIR"
    exit 1
fi

echo "Collecting best-model .cif file for each candidate..."
echo "=============================================="

# For each candidate subfolder, pick the model with the highest
# ranking_score (same logic as structure_qc_en.py) and copy it
# with a clean, unique filename into a flat folder Foldseek can
# scan directly.
python3 - "$STRUCTURE_DIR" "$DB_DIR/structures" << 'PYEOF'
import sys, os, glob, json, shutil

structure_dir, out_dir = sys.argv[1], sys.argv[2]

subfolders = [f for f in glob.glob(os.path.join(structure_dir, "*")) if os.path.isdir(f)]

for folder in sorted(subfolders):
    name = os.path.basename(folder)
    summary_files = glob.glob(os.path.join(folder, "*summary_confidences_*.json"))

    best_index, best_score = 0, float("-inf")
    for path in summary_files:
        base = os.path.basename(path)
        try:
            idx = int(base.rsplit("_", 1)[-1].replace(".json", ""))
        except ValueError:
            continue
        try:
            with open(path) as f:
                score = json.load(f).get("ranking_score", 0)
        except Exception:
            score = 0
        if score > best_score:
            best_score, best_index = score, idx

    cif_matches = glob.glob(os.path.join(folder, f"*model_{best_index}.cif"))
    if not cif_matches:
        cif_matches = glob.glob(os.path.join(folder, "*.cif"))
    if not cif_matches:
        print(f"  [!] {name}: no .cif file found, skipping.")
        continue

    dest = os.path.join(out_dir, f"{name}.cif")
    shutil.copy(cif_matches[0], dest)
    print(f"  {name}: using model_{best_index}.cif (ranking_score={best_score:.3f})")
PYEOF

N_STRUCTURES=$(ls "$DB_DIR/structures"/*.cif 2>/dev/null | wc -l)
echo ""
echo "Collected $N_STRUCTURES structures for comparison."
echo "=============================================="

if [ "$N_STRUCTURES" -eq 0 ]; then
    echo "ERROR: no structures collected, nothing to compare."
    exit 1
fi

echo "Building Foldseek database..."
foldseek createdb "$DB_DIR/structures" "$DB_DIR/db"

echo "Running all-vs-all structural comparison..."
mkdir -p "$OUT_DIR/tmp"
# NOTE: the -a flag is required to compute the full alignment
# backtrace, which is needed later to report TM-score columns
# in convertalis. Without it, convertalis fails/crashes when
# asked for qtmscore/ttmscore/alntmscore.
foldseek search "$DB_DIR/db" "$DB_DIR/db" "$OUT_DIR/results" "$OUT_DIR/tmp" \
    --alignment-type 1 -e 0.001 -a

echo "Converting results to a readable table..."
foldseek convertalis "$DB_DIR/db" "$DB_DIR/db" "$OUT_DIR/results" \
    "$OUT_DIR/structure_comparison.tsv" \
    --format-output "query,target,fident,alnlen,evalue,bits,qtmscore,ttmscore,alntmscore"

echo "=============================================="
echo "Foldseek comparison finished."
echo "Results in: $OUT_DIR/structure_comparison.tsv"
echo "=============================================="
echo ""
echo "Column meaning:"
echo "  fident      : sequence identity (fraction, 0-1)"
echo "  qtmscore/ttmscore/alntmscore : TM-score (structural similarity, 0-1)"
echo "                  >0.5 = same fold, >0.7 = highly similar structure"
