#!/bin/bash
# ==========================================================
# Combines qc_summary.csv (Phase 1) + nlr_summary.csv (Phase 2)
# into a single table, matching by genome name, and adds
# columns normalized by genome size:
#
#   - genome_size_mb           : genome size in Mb
#   - nlr_total_por_mb         : total NLR candidates / Mb
#   - nlr_completos_por_mb     : complete NLR (CNL+TNL) / Mb
#
# Normalizing by size matters because comparing the raw number
# of NLR genes between genomes of very different sizes
# (e.g. 328 Mb vs 2.7 Gb) would not be a fair comparison.
#
# Works even if Phase 2 (NLR) hasn't finished yet -- it simply
# combines whichever genomes are already present in both
# tables.
#
# Usage:
#   cd ~/project
#   bash combine_qc_nlr.sh
# ==========================================================

set -uo pipefail

QC_CSV="01_qc/qc_summary.csv"
NLR_CSV="02_nlr_annotation/nlr_summary.csv"
OUT_CSV="qc_nlr_combined.csv"

if [ ! -f "$QC_CSV" ]; then
    echo "ERROR: cannot find $QC_CSV"
    exit 1
fi

if [ ! -f "$NLR_CSV" ]; then
    echo "ERROR: cannot find $NLR_CSV"
    echo "Run first: cd 02_nlr_annotation && bash nlr_summary.sh"
    exit 1
fi

awk -F',' '
    NR==FNR {
        if (FNR==1) {
            qc_header = $0
            next
        }
        qc_data[$1] = $0
        qc_sumlen[$1] = $3
        next
    }
    FNR==1 {
        nlr_header = $0
        print qc_header "," nlr_header ",genome_size_mb,nlr_total_per_mb,nlr_complete_per_mb"
        next
    }
    {
        genome = $1
        gsub("/", "", genome)
        if (genome in qc_data) {
            sumlen = qc_sumlen[genome]
            size_mb = sumlen / 1000000
            total_candidates = $2
            nlr_complete = $5
            total_per_mb = total_candidates / size_mb
            complete_per_mb = nlr_complete / size_mb
            printf "%s,%s,%.2f,%.3f,%.3f\n", qc_data[genome], $0, size_mb, total_per_mb, complete_per_mb
        }
    }
' "$QC_CSV" FS=',' "$NLR_CSV" > "$OUT_CSV"

NUM_ROWS=$(( $(wc -l < "$OUT_CSV") - 1 ))
echo "Combined file saved to: $(pwd)/$OUT_CSV"
echo "Genomes combined: $NUM_ROWS"
echo ""
echo "New columns added:"
echo "  - genome_size_mb        : genome size in Mb"
echo "  - nlr_total_per_mb      : total NLR candidates / Mb of genome"
echo "  - nlr_complete_per_mb   : complete NLR (CC/TIR-NBARC-LRR) / Mb"
