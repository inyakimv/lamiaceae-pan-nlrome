#!/bin/bash
# ==========================================================
# NLR candidate summary per genome -> CSV for Excel
#
# Walks through all subfolders of 02_nlr_annotation/ and counts,
# for each genome, how many NLR candidates there are of each
# domain architecture (CC-NBARC-LRR, NBARC-LRR, TIR-NBARC-LRR,
# etc.)
#
# Usage:
#   cd ~/project/02_nlr_annotation
#   bash nlr_summary.sh
#
# Can be re-run at any time, even if not all 46 genomes are
# finished yet -- it only summarizes the ones already done.
# ==========================================================

set -uo pipefail

OUT_CSV="nlr_summary.csv"

echo "genome,total_candidates,CC_NBARC_LRR,TIR_NBARC_LRR,complete_NLR,CC_NBARC,NBARC_LRR,NBARC,TIR_NBARC,other" > "$OUT_CSV"

for dir in */; do
    genome="${dir%/}"
    nlr_file=$(find "$dir" -name "*_nlr.txt" 2>/dev/null | head -n1)

    if [ -z "$nlr_file" ] || [ ! -s "$nlr_file" ]; then
        continue
    fi

    total=$(wc -l < "$nlr_file")
    cc_nbarc_lrr=$(cut -f3 "$nlr_file" | grep -c "^CC-NBARC-LRR$")
    tir_nbarc_lrr=$(cut -f3 "$nlr_file" | grep -c "^TIR-NBARC-LRR$")
    cc_nbarc=$(cut -f3 "$nlr_file" | grep -c "^CC-NBARC$")
    nbarc_lrr=$(cut -f3 "$nlr_file" | grep -c "^NBARC-LRR$")
    nbarc=$(cut -f3 "$nlr_file" | grep -c "^NBARC$")
    tir_nbarc=$(cut -f3 "$nlr_file" | grep -c "^TIR-NBARC$")

    complete=$((cc_nbarc_lrr + tir_nbarc_lrr))
    known=$((cc_nbarc_lrr + tir_nbarc_lrr + cc_nbarc + nbarc_lrr + nbarc + tir_nbarc))
    other=$((total - known))

    echo "$genome,$total,$cc_nbarc_lrr,$tir_nbarc_lrr,$complete,$cc_nbarc,$nbarc_lrr,$nbarc,$tir_nbarc,$other" >> "$OUT_CSV"
done

echo "Summary saved to: $(pwd)/$OUT_CSV"
column -s',' -t "$OUT_CSV"
