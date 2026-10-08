#!/bin/bash
# ==========================================================
# Runs OrthoFinder on the 46 proteomes to group NLR genes
# into homologous gene families (orthogroups).
#
# IMPORTANT: run from ~/project
#
# Requirements (install once, in its OWN conda environment --
# installing into "base" makes the solver take forever
# because it has to check compatibility with everything else
# already installed):
#   conda create -n orthofinder_env --solver=classic -c bioconda -c conda-forge orthofinder -y
#
# Usage:
#   cd ~/project
#   bash run_orthofinder.sh
# ==========================================================

set -uo pipefail

PROTEOMES_DIR="03_orthofinder/proteomes"
THREADS=32
CONDA_ENV="orthofinder_env"

if [ ! -d "$PROTEOMES_DIR" ]; then
    echo "ERROR: cannot find folder $PROTEOMES_DIR"
    echo "Make sure you're running this from ~/project"
    exit 1
fi

N_FASTA=$(ls "$PROTEOMES_DIR"/*.fasta 2>/dev/null | wc -l)
if [ "$N_FASTA" -eq 0 ]; then
    echo "ERROR: no .fasta files found in $PROTEOMES_DIR"
    exit 1
fi

CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

echo "Found $N_FASTA proteomes in $PROTEOMES_DIR"
echo "Running OrthoFinder with $THREADS threads..."
echo "=============================================="

orthofinder -f "$PROTEOMES_DIR" -t "$THREADS" -a "$THREADS"

conda deactivate

echo "=============================================="
echo "OrthoFinder finished."
echo "Results in: $PROTEOMES_DIR/OrthoFinder/Results_<date>/"
echo "=============================================="
