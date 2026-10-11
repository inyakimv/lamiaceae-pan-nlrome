#!/bin/bash
# ==========================================================
# NLR domain search (HMMER) across ANNEVO-predicted proteomes
# for the 46 Lamiaceae species, using a reduced Pfam profile
# set (NB-ARC, LRR, TIR).
#
# Usage:
#   cd ~/project/annevo_proteins
#   bash run_hmmer_nlr_search.sh
# ==========================================================
set -uo pipefail

PROFILES="../nlr_domain_search/nlr_profiles.hmm"
OUT_DIR="../nlr_domain_search/per_species"
THREADS=8

mkdir -p "$OUT_DIR"

shopt -s nullglob
FASTA_FILES=(*.annevo_proteins.fasta)
if [ ${#FASTA_FILES[@]} -eq 0 ]; then
    echo "ERROR: no *.annevo_proteins.fasta files found in this folder."
    exit 1
fi

echo "Found ${#FASTA_FILES[@]} proteomes."

for FASTA in "${FASTA_FILES[@]}"; do
    BASENAME="${FASTA%.annevo_proteins.fasta}"
    OUT_DOMTBL="$OUT_DIR/${BASENAME}.domtblout"

    echo "=============================================="
    echo ">> Processing $BASENAME"
    echo "=============================================="

    if [ -f "$OUT_DOMTBL" ] && [ -s "$OUT_DOMTBL" ]; then
        echo "  -> $OUT_DOMTBL already exists, skipping."
        continue
    fi

    hmmscan --cpu "$THREADS" --domtblout "$OUT_DOMTBL" \
        "$PROFILES" "$FASTA" > /dev/null

    if [ -s "$OUT_DOMTBL" ]; then
        N_HITS=$(grep -vc "^#" "$OUT_DOMTBL")
        N_NBARC=$(grep -v "^#" "$OUT_DOMTBL" | awk '$1=="NB-ARC"' | awk '{print $4}' | sort -u | wc -l)
        echo "  -> Done: $BASENAME ($N_HITS total domain hits, $N_NBARC proteins with NB-ARC)"
    else
        echo "  [!] No output generated for $BASENAME."
    fi
    echo ""
done

echo "=============================================="
echo "HMMER NLR domain search finished."
echo "Results in: $OUT_DIR/"
echo "=============================================="
