#!/bin/bash
# ==========================================================
# BLASTs the protein sequences of the highlighted orthogroups
# against the RefPlantNLR database, to find their closest
# experimentally validated NLR homolog.
#
# Up to 5 target sequences are reported per query; the summary
# script keeps, for every query, the alignment with the highest
# bit score. (-max_target_seqs 1 can return a hit that is not
# the best one.)
#
# IMPORTANT: run from ~/project
#
# Requirements: BLAST+ (tested with 2.12.0+)
#
# Usage:
#   cd ~/project
#   bash scripts/validation/run_blast_validation_en.sh OG0000000 OG0000008 OG0000018 OG0000029 OG0000038
# ==========================================================

set -uo pipefail

REF_DB_FASTA="refplantnlr/RefPlantNLR_proteins.fasta"
BLAST_DB="refplantnlr/RefPlantNLR_blastdb"
OUT_DIR="10_validation"
EVALUE="1e-10"

if [ ! -f "$REF_DB_FASTA" ]; then
    echo "ERROR: cannot find $REF_DB_FASTA"
    exit 1
fi

mkdir -p "$OUT_DIR"
blastp -version | head -1

# Build the BLAST database only once (skip if already built)
if [ ! -f "${BLAST_DB}.pin" ]; then
    echo "Building BLAST database from RefPlantNLR..."
    makeblastdb -in "$REF_DB_FASTA" -dbtype prot -out "$BLAST_DB"
fi

for og_id in "$@"; do
    query="${OUT_DIR}/${og_id}_proteins_for_blast.fasta"
    output="${OUT_DIR}/${og_id}_blast_results.tsv"

    if [ ! -f "$query" ]; then
        echo "[!] Cannot find $query for $og_id, skipping."
        echo "    Run first: python3 scripts/structure/extract_proteins_for_blast_en.py $og_id"
        continue
    fi

    echo "Running BLAST for $og_id..."
    blastp -query "$query" -db "$BLAST_DB" \
        -evalue "$EVALUE" \
        -outfmt "6 qseqid sseqid pident length qlen slen evalue bitscore stitle" \
        -max_target_seqs 5 \
        -num_threads 8 \
        -out "$output"

    n_hits=$(cut -f1 "$output" | sort -u | wc -l)
    n_total=$(grep -c "^>" "$query")
    echo "  -> $og_id: $n_hits / $n_total sequences with a hit (e-value <= $EVALUE)"
done

echo "=============================================="
echo "BLAST finished. Results in: $OUT_DIR/<OG>_blast_results.tsv"
echo "Next: python3 scripts/validation/summarize_blast_results_en.py"
echo "=============================================="
