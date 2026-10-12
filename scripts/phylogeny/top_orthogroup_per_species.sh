#!/bin/bash
# ==========================================================
# For each of the 46 species, finds its most expanded
# orthogroup and calculates what percentage of that group's
# total copies belong to it.
#
# Useful for automatically identifying, across all 46
# species, which ones have expansions concentrated in few
# groups (candidates for deeper tandem duplication analysis).
#
# IMPORTANT: run from the Orthogroups/ folder
#
# Usage:
#   cd ~/project/03_orthofinder/proteomes/OrthoFinder/Results_Jul07/Orthogroups
#   bash top_orthogroup_per_species.sh
# ==========================================================

set -uo pipefail

INPUT="Orthogroups.GeneCount.tsv"
OUTPUT="top_orthogroup_summary.csv"

echo "species,top_orthogroup,copies_in_species,total_in_group,pct_concentration" > "$OUTPUT"

# Read header to know species column names/positions
HEADER=$(head -1 "$INPUT")
NCOLS=$(echo "$HEADER" | awk -F'\t' '{print NF}')

# Loop over each species column (from 2 to second-to-last,
# since column 1 is "Orthogroup" and the last is "Total")
for col in $(seq 2 $((NCOLS-1))); do
    species=$(echo "$HEADER" | awk -F'\t' -v c=$col '{print $c}')

    awk -F'\t' -v c=$col -v sp="$species" '
        NR==1 {next}
        {
            if ($c > max) {
                max = $c
                best_og = $1
                best_total = $NF
            }
        }
        END {
            pct = (best_total > 0) ? (max/best_total*100) : 0
            printf "%s,%s,%d,%d,%.1f\n", sp, best_og, max, best_total, pct
        }
    ' "$INPUT" >> "$OUTPUT"
done

echo "Summary saved to: $(pwd)/$OUTPUT"
column -s',' -t "$OUTPUT" | sort -t' ' -k5 -rn | head -50
