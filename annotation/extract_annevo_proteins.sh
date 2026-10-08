#!/bin/bash
# ==========================================================
# Protein extraction from ANNEVO gene annotations (46 species)
# using AGAT (agat_sp_extract_sequences.pl).
#
# Requirements:
#   conda activate ANNEVO   (AGAT is installed inside this env)
#
# Usage:
#   cd ~/project/00_genomes
#   bash extract_annevo_proteins.sh
# ==========================================================
set -uo pipefail

OUT_DIR="../annevo_proteins"
mkdir -p "$OUT_DIR"

shopt -s nullglob
GFF_FILES=(*.annevo.gff3)
if [ ${#GFF_FILES[@]} -eq 0 ]; then
    echo "ERROR: no *.annevo.gff3 files found in this folder."
    exit 1
fi

echo "Found ${#GFF_FILES[@]} ANNEVO GFF3 files."

for GFF in "${GFF_FILES[@]}"; do
    BASENAME="${GFF%.annevo.gff3}"
    FNA="${BASENAME}.fna"
    FNA_GZ="${BASENAME}.fna.gz"
    OUT_PROT="$OUT_DIR/${BASENAME}.annevo_proteins.fasta"

    echo "=============================================="
    echo ">> Processing $BASENAME"
    echo "=============================================="

    # Skip if a valid protein file already exists
    if [ -f "$OUT_PROT" ] && [ -s "$OUT_PROT" ]; then
        N_PROT=$(grep -c "^>" "$OUT_PROT" 2>/dev/null || echo 0)
        if [ "$N_PROT" -gt 0 ]; then
            echo "  -> $OUT_PROT already exists with $N_PROT proteins, skipping."
            continue
        fi
    fi

    # Check available disk space before each genome
    AVAIL_KB=$(df --output=avail "$HOME" | tail -1)
    AVAIL_GB=$((AVAIL_KB / 1024 / 1024))
    if [ "$AVAIL_GB" -lt 10 ]; then
        echo "  [!] WARNING: less than 10G free on disk (${AVAIL_GB}G). Aborting."
        break
    fi

    if [ ! -f "$FNA" ]; then
        if [ ! -f "$FNA_GZ" ]; then
            echo "  [!] Neither $FNA nor $FNA_GZ found. Skipping."
            continue
        fi
        gunzip -k "$FNA_GZ"
    fi

    agat_sp_extract_sequences.pl -g "$GFF" -f "$FNA" -p -o "$OUT_PROT" 2>&1 | tail -20

    if [ -s "$OUT_PROT" ]; then
        N_PROT=$(grep -c "^>" "$OUT_PROT" 2>/dev/null || echo 0)
        echo "  -> Done: $BASENAME ($N_PROT proteins extracted)"
    else
        echo "  [!] No protein output generated for $BASENAME."
    fi

    # Cleanup: remove decompressed .fna (regenerable from .gz)
    rm -f "$FNA"
    echo ""
done

echo "=============================================="
echo "Protein extraction finished."
echo "Results in: $OUT_DIR/"
echo "=============================================="
