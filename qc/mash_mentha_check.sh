#!/bin/bash
# Genome-level sanity check of the Mentha placement in the species tree.
# Tool: Mash v2.3 (Ondov et al., 2016), conda env "mash_env":
#   conda create -n mash_env -c bioconda -c conda-forge mash -y
# Output: results_tables/mash_mentha_check.tsv  (Mash distances, all pairs)
set -euo pipefail
PROJECT="${PROJECT:-$HOME/project}"
cd "$PROJECT/00_genomes"
mkdir -p "$PROJECT/results_tables"
conda run -n mash_env mash sketch -p 16 -k 21 -s 1000000 -o "$PROJECT/mentha_check" \
  *Mentha*_genomic.fna.gz *Drepanocaryum*_genomic.fna.gz *Lycopus*_genomic.fna.gz *Salvia_rosmarinus*_genomic.fna.gz
conda run -n mash_env mash dist "$PROJECT/mentha_check.msh" "$PROJECT/mentha_check.msh" \
  > "$PROJECT/results_tables/mash_mentha_check.tsv"
