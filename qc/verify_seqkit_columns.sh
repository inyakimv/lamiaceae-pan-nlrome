#!/bin/bash
cd ~/project/00_genomes
OUT=~/project/seqkit_check.csv
touch "$OUT"
for f in *_genomic.fna.gz; do
  b=${f%.fna.gz}
  grep -q "^$b," "$OUT" && continue
  seqkit stats -a -T "$f" | tail -n1 | awk -F'\t' -v b="$b" '{print b","$4","$5","$13","$18}' >> "$OUT"
done
diff <(sort "$OUT") <(tail -n +2 ~/project/01_qc/qc_summary.csv | cut -d, -f1-5 | sort) && echo IDENTICAL
