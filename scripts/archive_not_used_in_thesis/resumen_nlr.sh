#!/bin/bash
# ==========================================================
# Resumen de candidatos NLR por genoma -> CSV para Excel
#
# Recorre todas las subcarpetas de 02_nlr_annotation/ y cuenta,
# por cada genoma, cuantos candidatos NLR hay de cada tipo
# (CC-NBARC-LRR, NBARC-LRR, TIR-NBARC-LRR, etc.)
#
# Uso:
#   cd ~/project/02_nlr_annotation
#   bash resumen_nlr.sh
#
# Se puede volver a correr en cualquier momento, aunque no
# hayan terminado los 46 genomas -- solo resume los que ya
# esten listos.
# ==========================================================

set -uo pipefail

OUT_CSV="nlr_summary.csv"

echo "genoma,total_candidatos,CC_NBARC_LRR,TIR_NBARC_LRR,NLR_completos,CC_NBARC,NBARC_LRR,NBARC,TIR_NBARC,otros" > "$OUT_CSV"

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

    completos=$((cc_nbarc_lrr + tir_nbarc_lrr))
    conocidos=$((cc_nbarc_lrr + tir_nbarc_lrr + cc_nbarc + nbarc_lrr + nbarc + tir_nbarc))
    otros=$((total - conocidos))

    echo "$genome,$total,$cc_nbarc_lrr,$tir_nbarc_lrr,$completos,$cc_nbarc,$nbarc_lrr,$nbarc,$tir_nbarc,$otros" >> "$OUT_CSV"
done

echo "Resumen guardado en: $(pwd)/$OUT_CSV"
column -s',' -t "$OUT_CSV"
