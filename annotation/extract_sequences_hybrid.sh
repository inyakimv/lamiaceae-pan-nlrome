#!/bin/bash
# ==========================================================
# Extracts protein sequences of "complete" NLR genes using a
# HYBRID approach:
#   - Primary method: miniprot, aligning against RefPlantNLR
#   - Fallback method: getorf (direct DNA translation) for any
#     locus where miniprot found NO alignment
#
# v5: fixes a bug where getorf mangles FASTA headers that
# contain "::" (EMBOSS treats it as a database:id separator
# and silently discards everything before the last colon).
# Headers are now cleaned to plain locus IDs before being
# passed to getorf.
#
# IMPORTANT: run from ~/project
#
# Requirements (install once):
#   conda install --solver=classic -c bioconda bedtools miniprot samtools emboss
#   pip install biopython --break-system-packages
#
# Usage:
#   cd ~/project
#   bash extract_sequences_hybrid_v5.sh
# ==========================================================

set -uo pipefail

GENOMES_DIR="00_genomes"
NLR_DIR="02_nlr_annotation"
OUT_DIR="02c_nlr_refined"
ORTHO_DIR="03_orthofinder/proteomes"
REF_DB="refplantnlr/RefPlantNLR_proteins.fasta"
FLANK_BP=1500
THREADS=8

mkdir -p "$OUT_DIR"
mkdir -p "$ORTHO_DIR"

if [ ! -f "$REF_DB" ]; then
    echo "ERROR: cannot find reference database $REF_DB"
    exit 1
fi

for dir in "$NLR_DIR"/*/; do
    genome=$(basename "$dir")
    nlr_txt=$(find "$dir" -name "*_nlr.txt" 2>/dev/null | head -n1)
    fna="$GENOMES_DIR/${genome}.fna"

    if [ -z "$nlr_txt" ] || [ ! -s "$nlr_txt" ]; then
        echo "  [!] No NLR results for $genome. Skipping."
        continue
    fi
    if [ ! -f "$fna" ]; then
        echo "  [!] Cannot find genome $fna. Skipping."
        continue
    fi

    echo "=============================================="
    echo ">> Processing $genome"
    echo "=============================================="

    WORKDIR="$OUT_DIR/$genome"
    mkdir -p "$WORKDIR"

    if [ ! -f "${fna}.fai" ]; then
        samtools faidx "$fna"
    fi

    awk -F'\t' '$3=="CC-NBARC-LRR" || $3=="TIR-NBARC-LRR" {
        printf "%s\t%d\t%s\t%s\t0\t%s\n", $1, $4-1, $5, $2, $6
    }' "$nlr_txt" > "$WORKDIR/nlr_complete.bed"

    N_LOCI=$(wc -l < "$WORKDIR/nlr_complete.bed")
    if [ "$N_LOCI" -eq 0 ]; then
        echo "  [!] $genome has no complete NLR. Skipping."
        continue
    fi
    echo "  -> $N_LOCI complete NLR loci"

    # --- Track A: miniprot (flanked region, aligned to RefPlantNLR) ---
    bedtools slop -i "$WORKDIR/nlr_complete.bed" -g "${fna}.fai" \
        -b "$FLANK_BP" > "$WORKDIR/nlr_flanked.bed"

    bedtools getfasta -fi "$fna" -bed "$WORKDIR/nlr_flanked.bed" \
        -s -name -fo "$WORKDIR/nlr_flanked.fasta" 2> "$WORKDIR/bedtools_flanked_log.txt"

    miniprot --gff -t "$THREADS" "$WORKDIR/nlr_flanked.fasta" "$REF_DB" \
        > "$WORKDIR/miniprot.gff" 2> "$WORKDIR/miniprot_log.txt"

    # --- Track B: getorf fallback (plain, non-flanked region) ---
    bedtools getfasta -fi "$fna" -bed "$WORKDIR/nlr_complete.bed" \
        -s -name -fo "$WORKDIR/nlr_plain.fasta" 2> "$WORKDIR/bedtools_plain_log.txt"

    # IMPORTANT FIX: clean headers to plain locus_id before getorf,
    # since getorf/EMBOSS mangles headers containing "::"
    awk '/^>/{split($0,a,"::"); print a[1]; next} {print}' \
        "$WORKDIR/nlr_plain.fasta" > "$WORKDIR/nlr_plain_clean.fasta"

    getorf -sequence "$WORKDIR/nlr_plain_clean.fasta" \
        -outseq "$WORKDIR/nlr_plain_orfs.fasta" \
        -find 1 -minsize 300 2> "$WORKDIR/getorf_log.txt"

    # --- Combine both tracks, preferring miniprot, falling back to getorf ---
    python3 - "$WORKDIR/miniprot.gff" "$WORKDIR/nlr_flanked.fasta" \
        "$WORKDIR/nlr_plain_orfs.fasta" "$WORKDIR/nlr_complete.bed" \
        "$WORKDIR/nlr_complete_proteins.fasta" << 'PYEOF'
import sys
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

gff_file, flanked_fasta, orfs_fasta, bed_file, out_file = sys.argv[1:6]

all_loci = []
with open(bed_file) as f:
    for line in f:
        cols = line.rstrip("\n").split("\t")
        if len(cols) >= 4:
            all_loci.append(cols[3])

# ---- Track A: parse miniprot results ----
seqs = SeqIO.to_dict(SeqIO.parse(flanked_fasta, "fasta"))
mrna_strand = {}
cds_by_mrna = defaultdict(list)

with open(gff_file) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        cols = line.rstrip("\n").split("\t")
        if len(cols) < 9:
            continue
        seqid, source, feature, start, end, score, strand, phase, attrs = cols
        if feature == "mRNA":
            attr_dict = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
            mrna_id = attr_dict.get("ID", "")
            mrna_strand[mrna_id] = (seqid, strand)
        elif feature == "CDS":
            attr_dict = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
            parent = attr_dict.get("Parent", "")
            ph = int(phase) if phase.isdigit() else 0
            cds_by_mrna[parent].append((int(start), int(end), ph))

candidates_by_locus = defaultdict(list)

for mrna_id, cds_list in cds_by_mrna.items():
    if mrna_id not in mrna_strand:
        continue
    seqid, strand = mrna_strand[mrna_id]
    if seqid not in seqs:
        continue

    locus_id = seqid.split("::")[0]

    cds_list.sort(key=lambda x: x[0])
    full_seq = seqs[seqid].seq
    first_phase = cds_list[-1][2] if strand == "-" else cds_list[0][2]

    concat = "".join(str(full_seq[start-1:end]) for start, end, ph in cds_list)
    seq_obj = Seq(concat)
    if strand == "-":
        seq_obj = seq_obj.reverse_complement()
    if first_phase > 0:
        seq_obj = seq_obj[first_phase:]
    trim = len(seq_obj) % 3
    if trim:
        seq_obj = seq_obj[:-trim]

    protein = str(seq_obj.translate(to_stop=True))
    candidates_by_locus[locus_id].append(protein)

miniprot_proteins = {}
for locus_id, prots in candidates_by_locus.items():
    best = max(prots, key=len)
    if len(best) >= 100:
        miniprot_proteins[locus_id] = best

# ---- Track B: parse getorf fallback results (clean headers, no "::") ----
orf_seqs = defaultdict(dict)
current_id = None
current_seq = []

def flush():
    if current_id is not None:
        locus_id = current_id.rsplit("_", 1)[0]
        orf_seqs[locus_id][current_id] = "".join(current_seq)

with open(orfs_fasta) as f:
    for line in f:
        line = line.rstrip()
        if line.startswith(">"):
            flush()
            current_id = line[1:].split()[0]
            current_seq = []
        else:
            current_seq.append(line)
    flush()

getorf_proteins = {}
for locus_id, orfs in orf_seqs.items():
    longest_id = max(orfs, key=lambda k: len(orfs[k]))
    protein = orfs[longest_id]
    if len(protein) >= 100:
        getorf_proteins[locus_id] = protein

# ---- Combine: prefer miniprot, fall back to getorf ----
results = []
n_miniprot = 0
n_getorf = 0
n_none = 0

for locus_id in all_loci:
    if locus_id in miniprot_proteins:
        results.append((locus_id, miniprot_proteins[locus_id], "miniprot"))
        n_miniprot += 1
    elif locus_id in getorf_proteins:
        results.append((locus_id, getorf_proteins[locus_id], "getorf"))
        n_getorf += 1
    else:
        n_none += 1

with open(out_file, "w") as out:
    for locus_id, protein, source in results:
        out.write(f">{locus_id}|{source}\n{protein}\n")

print(f"  -> {len(results)} proteins total ({n_miniprot} via miniprot, {n_getorf} via getorf fallback, {n_none} could not be translated by either method)")
PYEOF

    N_PROT=$(grep -c "^>" "$WORKDIR/nlr_complete_proteins.fasta" 2>/dev/null || echo 0)
    echo "  -> Total: $N_PROT protein sequences (hybrid)"

    short_name=$(echo "$genome" | cut -d'_' -f1,2)
    sed "s/^>/>${short_name}|/" "$WORKDIR/nlr_complete_proteins.fasta" \
        > "$ORTHO_DIR/${genome}.fasta"

    echo ""
done

echo "=============================================="
echo "Hybrid extraction finished."
echo "Proteomes ready for OrthoFinder in: $ORTHO_DIR/"
echo "=============================================="
