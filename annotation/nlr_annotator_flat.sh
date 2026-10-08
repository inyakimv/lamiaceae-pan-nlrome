#!/bin/bash
# ==========================================================
# NLR gene annotation with NLR-Annotator
# Version adapted for a FLAT folder (all .fna.gz files
# directly in the current directory, no subfolders)
#
# IMPORTANT: this script must be run from a folder that is a
# DIRECT CHILD of ~/project (e.g. ~/project/00_genomes or
# ~/project/nlr_test), because it saves results to
# ../02_nlr_annotation/ (one level up)
#
# Requirements (install once, see STEP 2.1 of the guide):
#   conda install -c conda-forge openjdk
#   git clone -b NLR-Annotator-2 https://github.com/steuernb/NLR-Annotator.git ~/software/NLR-Annotator
#
# Usage:
#   cd ~/project/00_genomes   (or ~/project/nlr_test)
#   bash nlr_annotator_flat.sh
# ==========================================================

set -uo pipefail

# Adjust these paths if you installed NLR-Annotator elsewhere
NLR_JAR="$HOME/software/NLR-Annotator/NLR-Annotator-v2.1b.jar"
MOT_FILE="$HOME/software/NLR-Annotator/src/mot.txt"
STORE_FILE="$HOME/software/NLR-Annotator/src/store.txt"
THREADS=8
JAVA_MEM="8000M"

OUT_DIR="../02_nlr_annotation"
mkdir -p "$OUT_DIR"

if [ ! -d "$OUT_DIR" ]; then
    echo "ERROR: cannot find folder $OUT_DIR"
    echo "Make sure you're running this script from a folder that"
    echo "is a direct child of ~/project (e.g: 00_genomes or nlr_test)."
    exit 1
fi

for f in "$NLR_JAR" "$MOT_FILE" "$STORE_FILE"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: cannot find file: $f"
        echo "Check STEP 2.1 of the guide (NLR-Annotator installation)"
        echo "and adjust the NLR_JAR / MOT_FILE / STORE_FILE variables in this script."
        exit 1
    fi
done

shopt -s nullglob
FILES=(*_genomic.fna.gz)

if [ ${#FILES[@]} -eq 0 ]; then
    echo "ERROR: no *_genomic.fna.gz files found in this folder."
    exit 1
fi

echo "Found ${#FILES[@]} genomes."
echo "Results will be saved to: $OUT_DIR"

for FNA_GZ in "${FILES[@]}"; do
    BASENAME="${FNA_GZ%.fna.gz}"
    echo "=============================================="
    echo ">> Processing $BASENAME"
    echo "=============================================="

    FNA="${BASENAME}.fna"
    if [ ! -f "$FNA" ]; then
        gunzip -k "$FNA_GZ"
    fi

    GENOME_OUT="$OUT_DIR/$BASENAME"
    mkdir -p "$GENOME_OUT"

    java -Xmx${JAVA_MEM} -jar "$NLR_JAR" \
        -i "$FNA" \
        -x "$MOT_FILE" \
        -y "$STORE_FILE" \
        -o "$GENOME_OUT/${BASENAME}_nlr.txt" \
        -g "$GENOME_OUT/${BASENAME}_nlr.gff" \
        -b "$GENOME_OUT/${BASENAME}_nlr.bed" \
        -t "$THREADS" \
        > "$GENOME_OUT/${BASENAME}_nlr_log.txt" 2>&1

    if [ -s "$GENOME_OUT/${BASENAME}_nlr.txt" ]; then
        NUM_NLR=$(wc -l < "$GENOME_OUT/${BASENAME}_nlr.txt")
        echo "  -> Done: $BASENAME ($NUM_NLR NLR candidates found)"
    else
        echo "  [!] No output generated for $BASENAME. Check $GENOME_OUT/${BASENAME}_nlr_log.txt"
    fi
    echo ""
done

echo "=============================================="
echo "NLR annotation finished."
echo "Per-genome results in: $OUT_DIR/<genome_name>/"
echo "=============================================="
