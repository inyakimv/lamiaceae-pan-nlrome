#!/bin/bash
# ==========================================================
# Gene annotation with ANNEVO for the 46 Lamiaceae species
# GPU-based, automatic cleanup of temp files and tmp/tmp_*
# folders after each genome (ANNEVO does not clean these up
# on its own).
#
# Requirements:
#   conda activate ANNEVO
#   Run from inside a screen/tmux session (long-running job)
#
# Usage:
#   cd ~/project/00_genomes
#   bash run_annevo_batch.sh
# ==========================================================
set -uo pipefail

ANNEVO_SCRIPT="$HOME/software/01_Scripts/RunANNEVOv223.sh"
GPU_ID=1
THREADS=24
LINEAGE="Magnoliopsida"

mkdir -p ./tmp

shopt -s nullglob
FILES=(*_genomic.fna.gz)
if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: no *_genomic.fna.gz files found in this folder."
    exit 1
fi

echo "Found ${#FILES[@]} genomes."

for FNA_GZ in "${FILES[@]}"; do
    BASENAME="${FNA_GZ%.fna.gz}"
    FNA="${BASENAME}.fna"
    OUT_GFF="${BASENAME}.annevo.gff3"
    TEMP_GFF="temp_${OUT_GFF}"
    TEMP_FORMATTED="temp_formatted_${OUT_GFF}"

    echo "=============================================="
    echo ">> Processing $BASENAME"
    echo "=============================================="

    # Skip if a VALID result already exists (contains real genes)
    if [ -f "$OUT_GFF" ] && [ -s "$OUT_GFF" ]; then
        N_GENES=$(grep -c $'\tgene\t' "$OUT_GFF" 2>/dev/null || echo 0)
        if [ "$N_GENES" -gt 0 ]; then
            echo "  -> $OUT_GFF already exists with $N_GENES genes, skipping."
            continue
        else
            echo "  -> $OUT_GFF exists but is empty/invalid. Reprocessing."
        fi
    fi

    # Check available disk space before each genome (abort if <10G free)
    AVAIL_KB=$(df --output=avail "$HOME" | tail -1)
    AVAIL_GB=$((AVAIL_KB / 1024 / 1024))
    if [ "$AVAIL_GB" -lt 10 ]; then
        echo "  [!] WARNING: less than 10G free on disk (${AVAIL_GB}G). Aborting to avoid quota failures."
        break
    fi

    if [ ! -f "$FNA" ]; then
        gunzip -k "$FNA_GZ"
    fi

    "$ANNEVO_SCRIPT" -g "$FNA" -o "$OUT_GFF" -c "$GPU_ID" -t "$THREADS" -l "$LINEAGE"

    if [ -s "$OUT_GFF" ]; then
        N_GENES=$(grep -c $'\tgene\t' "$OUT_GFF" 2>/dev/null || echo 0)
        echo "  -> Done: $BASENAME ($N_GENES genes annotated)"
    else
        echo "  [!] No output generated for $BASENAME. Check the log."
    fi

    # Cleanup: temp gff3 files + decompressed .fna + tmp/tmp_* folders
    # (ANNEVO does not remove these on its own; they can accumulate
    # to hundreds of GB across 46 genomes if left untouched)
    rm -f "$TEMP_GFF" "$TEMP_FORMATTED"
    rm -f "$FNA"
    find ./tmp -maxdepth 1 -type d -name "tmp_*" -exec rm -rf {} +
    echo "  -> Temp files, decompressed .fna, and tmp/tmp_* folders removed."
    echo ""
done

echo "=============================================="
echo "ANNEVO annotation finished."
echo "=============================================="
