#!/bin/bash
# ==========================================================
# NB-ARC phylogeny of the Lamiaceae proteins of Section 4.4 together with
# RefPlantNLR (checks whether the BLAST best hit agrees with the position
# of each protein in a tree).
#
# IMPORTANT: run from ~/project
#
# Requirements: HMMER 3.3.2 (hmmscan), MAFFT, IQ-TREE 3 (iqtree3), Python 3
# + Biopython. Needs the profile database nlr_profiles.hmm (Pfam 37.1,
# already hmmpress-ed) built for Annex A7, and the files of Annex A20.
#
# Usage (three stages):
#   export NLR_HMM=/path/to/nlr_profiles.hmm        # your Pfam profile file
#   bash scripts/validation/run_nbarc_phylogeny_en.sh prepare   # fast
#   nohup bash scripts/validation/run_nbarc_phylogeny_en.sh tree > 11_nbarc_tree/tree.log 2>&1 &   # slow
#   bash scripts/validation/run_nbarc_phylogeny_en.sh placement # fast
# ==========================================================
set -euo pipefail

W="11_nbarc_tree"
QUERIES="10_validation/OG*_proteins_for_blast.fasta"
REF="refplantnlr/RefPlantNLR_proteins.fasta"
BLAST_SUMMARY="10_validation/blast_validation_summary.csv"
TOOLS="scripts/validation/nbarc_tools_en.py"
THREADS="${THREADS:-8}"
NLR_HMM="${NLR_HMM:-nlr_profiles.hmm}"
stage="${1:-}"

case "$stage" in
prepare)
    [ -f "$NLR_HMM" ] || { echo "ERROR: cannot find $NLR_HMM (set NLR_HMM=/path/to/nlr_profiles.hmm)"; exit 1; }
    [ -f "$REF" ] || { echo "ERROR: cannot find $REF"; exit 1; }
    mkdir -p "$W"
    { hmmscan -h | sed -n 2p; mafft --version 2>&1 | head -1; iqtree3 --version 2>&1 | head -1; } | tee "$W/versions.txt"

    python3 "$TOOLS" prepare --queries "$QUERIES" --ref "$REF" --outdir "$W"

    # NB-ARC domain of every protein: default thresholds of the Pfam
    # gathering cut-off (--cut_ga) of the NB-ARC profile
    hmmscan --cut_ga --cpu "$THREADS" --domtblout "$W/all_proteins.domtblout" \
        "$NLR_HMM" "$W/all_proteins.fasta" > /dev/null
    python3 "$TOOLS" extract --workdir "$W" --domtblout "$W/all_proteins.domtblout"

    mafft --auto --thread "$THREADS" "$W/nbarc_domains.fasta" > "$W/nbarc_aligned.fasta" 2> "$W/mafft.log"
    python3 "$TOOLS" trim --aln "$W/nbarc_aligned.fasta" --out "$W/nbarc_trimmed.fasta"
    echo "Next: nohup bash $0 tree > $W/tree.log 2>&1 &"
    ;;
tree)
    [ -f "$W/nbarc_trimmed.fasta" ] || { echo "ERROR: run the prepare stage first"; exit 1; }
    iqtree3 -s "$W/nbarc_trimmed.fasta" -m MFP -bb 1000 -nt 4 -pre "$W/nbarc_tree"
    echo "Tree: $W/nbarc_tree.treefile"
    ;;
placement)
    [ -f "$W/nbarc_tree.treefile" ] || { echo "ERROR: $W/nbarc_tree.treefile not found (tree stage not finished?)"; exit 1; }
    python3 "$TOOLS" placement --tree "$W/nbarc_tree.treefile" --idmap "$W/idmap.tsv" \
        --blast-summary "$BLAST_SUMMARY" --out "$W/nbarc_tree_placement.csv"
    ;;
*)
    echo "Usage: bash $0 prepare|tree|placement"; exit 1
    ;;
esac
