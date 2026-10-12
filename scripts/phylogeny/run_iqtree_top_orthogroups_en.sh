#!/bin/bash
# ==========================================================
# Runs IQ-TREE (with bootstrap) on the "top" orthogroups
# identified in Step 8/9 of Phase 3, reusing the alignments
# that OrthoFinder already generated.
#
# Instead of relying on OrthoFinder's fast trees (FastTree,
# without robust bootstrap), this rebuilds each tree with
# evolutionary model selection + 1000 ultrafast bootstrap
# (UFBoot) replicates -- the expected standard in
# phylogenetics for publication.
#
# IMPORTANT: run from the Results_<date>/ folder
#
# Requires: IQ-TREE (already installed as an OrthoFinder
# dependency, inside the orthofinder_env environment)
#
# Usage:
#   cd ~/project/03_orthofinder/proteomes/OrthoFinder/Results_Jul07
#   bash run_iqtree_top_orthogroups_en.sh
# ==========================================================

set -uo pipefail

MSA_DIR="MultipleSequenceAlignments"
OG_LIST_CSV="Orthogroups/top_orthogroup_summary.csv"
OUT_DIR="../../../../04_phylogeny"
THREADS=4
BOOTSTRAP=1000

# IMPORTANT: verify the actual name of the IQ-TREE executable
# in your installation before running this (could be iqtree,
# iqtree2, or iqtree3 depending on version). Adjust this
# variable:
IQTREE_BIN="iqtree3"

if ! command -v "$IQTREE_BIN" &> /dev/null; then
    echo "ERROR: cannot find command '$IQTREE_BIN'."
    echo "Check the correct name with: which iqtree iqtree2 iqtree3"
    echo "and adjust the IQTREE_BIN variable at the top of this script."
    exit 1
fi

if [ ! -f "$OG_LIST_CSV" ]; then
    echo "ERROR: cannot find $OG_LIST_CSV"
    echo "Run Step 8 of the OrthoFinder guide first."
    exit 1
fi

if [ ! -d "$MSA_DIR" ]; then
    echo "ERROR: cannot find folder $MSA_DIR"
    exit 1
fi

mkdir -p "$OUT_DIR"

# Extract the list of UNIQUE orthogroups (column 2 of the CSV,
# skipping header, removing duplicates)
UNIQUE_OGS=$(tail -n +2 "$OG_LIST_CSV" | cut -d',' -f2 | sort -u)
N_TOTAL=$(echo "$UNIQUE_OGS" | wc -l)

echo "Found $N_TOTAL unique orthogroups to process."
echo "Results will be saved to: $OUT_DIR"
echo "=============================================="

i=0
for og in $UNIQUE_OGS; do
    i=$((i+1))
    ALIGN_FILE="$MSA_DIR/${og}.fa"

    if [ ! -f "$ALIGN_FILE" ]; then
        echo "  [$i/$N_TOTAL] [!] Cannot find alignment for $og. Skipping."
        continue
    fi

    echo "  [$i/$N_TOTAL] Processing $og..."

    OG_OUT="$OUT_DIR/$og"
    mkdir -p "$OG_OUT"

    "$IQTREE_BIN" -s "$ALIGN_FILE" \
        -m MFP \
        -bb "$BOOTSTRAP" \
        -nt "$THREADS" \
        -pre "$OG_OUT/$og" \
        -quiet

    if [ -f "${OG_OUT}/${og}.treefile" ]; then
        echo "  [$i/$N_TOTAL] -> $og: tree generated successfully"
    else
        echo "  [$i/$N_TOTAL] [!] $og: tree not generated, check ${OG_OUT}/${og}.log"
    fi
done

echo "=============================================="
echo "IQ-TREE finished for $N_TOTAL orthogroups."
echo "Results in: $OUT_DIR/<orthogroup>/"
echo "=============================================="
