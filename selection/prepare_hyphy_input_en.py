#!/usr/bin/env python3
"""
Prepares the codon alignment (in-frame nucleotides) that HyPhy
needs, for ONE specific orthogroup.

For each protein in the group (already aligned by OrthoFinder),
recovers its corresponding nucleotide sequence:
  - If it came from miniprot: recomputes the CDS from the saved .gff
  - If it came from getorf: reruns getorf asking for nucleotides
    (-type N) and matches by length to the already-chosen protein

Then "threads" those nucleotides into the existing protein
alignment (respecting the gaps), producing a codon alignment
ready for HyPhy.

Usage:
    cd ~/project
    python3 prepare_hyphy_input_en.py OG0000038

Requires:
    - 03_orthofinder/proteomes/OrthoFinder/Results_*/MultipleSequenceAlignments/<OG>.fa
    - 02c_nlr_refined/<genome>/miniprot.gff and nlr_flanked.fasta (for miniprot)
    - 02c_nlr_refined/<genome>/nlr_plain_clean.fasta (for getorf)
    - EMBOSS (getorf) and Biopython installed
"""

import sys
import os
import glob
import subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

# Tolerance (in amino acids) between the protein length already
# chosen by OrthoFinder's input and the length of the recovered
# nucleotide ORF. Widened from 10 to 15 to rescue the closest
# candidates for species that would otherwise have ZERO
# representation in a given orthogroup (a representativeness
# issue, not just a "nice to have").
LENGTH_TOLERANCE_AA = 15

PROJECT_DIR = os.path.expanduser("~/project")
REFINED_DIR = os.path.join(PROJECT_DIR, "02c_nlr_refined")
ORTHOFINDER_RESULTS = glob.glob(os.path.join(
    PROJECT_DIR, "03_orthofinder/proteomes/OrthoFinder/Results_*"
))


def find_msa_file(og_id):
    for results_dir in ORTHOFINDER_RESULTS:
        candidate = os.path.join(results_dir, "MultipleSequenceAlignments", f"{og_id}.fa")
        if os.path.exists(candidate):
            return candidate
    return None


def parse_header(header):
    """Expected format: <something>|locus_id|method
    NOTE: OrthoFinder rewrites the part before the first '|' when
    it builds the alignment (it prefixes the sanitized filename,
    e.g. 'GCA_028654395_1_Scutellaria_barbata_genomic_GCA_028654395.1'
    instead of the plain 'GCA_028654395.1' we originally used). We
    extract the real accession (GCA_xxxxxxxxx.x / GCF_xxxxxxxxx.x)
    with a regex instead of assuming it's the whole first field."""
    import re
    parts = header.split("|")
    if len(parts) != 3:
        return None
    short_genome_field, locus_id, method = parts
    m = re.search(r'(GC[AF]_\d+\.\d+)', short_genome_field)
    if not m:
        return None
    return m.group(1), locus_id, method


def find_genome_folder(short_genome):
    """The header only carries a short prefix (e.g. GCF_023119035),
    but the real folder has the full genome name. We search by
    prefix match."""
    if not os.path.isdir(REFINED_DIR):
        return None
    for name in os.listdir(REFINED_DIR):
        if name.startswith(short_genome):
            return os.path.join(REFINED_DIR, name)
    return None


def get_nucleotide_from_miniprot(genome_folder, locus_id, target_protein_len):
    """Recomputes the nucleotide CDS from the saved miniprot gff,
    replicating the same phase logic used by the extraction script."""
    gff_file = os.path.join(genome_folder, "miniprot.gff")
    fasta_file = os.path.join(genome_folder, "nlr_flanked.fasta")
    if not os.path.exists(gff_file) or not os.path.exists(fasta_file):
        return None

    seqs = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))
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
            if seqid.split("::")[0] != locus_id:
                continue
            if feature == "mRNA":
                attr_dict = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
                mrna_id = attr_dict.get("ID", "")
                mrna_strand[mrna_id] = (seqid, strand)
            elif feature == "CDS":
                attr_dict = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
                parent = attr_dict.get("Parent", "")
                ph = int(phase) if phase.isdigit() else 0
                cds_by_mrna[parent].append((int(start), int(end), ph))

    for mrna_id, cds_list in cds_by_mrna.items():
        if mrna_id not in mrna_strand:
            continue
        seqid, strand = mrna_strand[mrna_id]
        if seqid not in seqs:
            continue

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

        # Only useful if it matches the length of the already-chosen
        # protein (recall the original script kept the best alignment
        # per locus; there may be more than one mRNA candidate per locus)
        translated_len = len(str(seq_obj.translate(to_stop=True)))
        diff = abs(translated_len - target_protein_len)
        if diff <= LENGTH_TOLERANCE_AA:
            nt = str(seq_obj)
            exact_len = target_protein_len * 3
            if len(nt) >= exact_len:
                return nt[:exact_len]
            return nt

    return None


def get_nucleotide_from_getorf(genome_folder, locus_id, target_protein_len):
    """Runs (or reuses) getorf in nucleotide mode on the clean
    input, and looks for the ORF whose amino acid length matches
    the already-chosen protein."""
    input_fasta = os.path.join(genome_folder, "nlr_plain_clean.fasta")
    if not os.path.exists(input_fasta):
        return None

    nucl_orfs_file = os.path.join(genome_folder, "nlr_plain_orfs_nucl.fasta")
    if not os.path.exists(nucl_orfs_file):
        # NOTE: this EMBOSS version (6.6.0.0) has no separate -type
        # flag; -find itself selects both output type and ORF
        # definition. -find 3 = nucleic sequences between START and
        # STOP codons (the nucleotide equivalent of the -find 1 used
        # for the protein extraction step).
        result = subprocess.run([
            "getorf", "-sequence", input_fasta,
            "-outseq", nucl_orfs_file,
            "-find", "3", "-minsize", "300"
        ], capture_output=True, text=True, check=False)
        if result.returncode != 0:
            print(f"  [!] getorf (nucleotide mode) failed for {genome_folder}: {result.stderr.strip()}")

    if not os.path.exists(nucl_orfs_file):
        return None

    candidates = []
    current_id, current_seq = None, []

    def flush():
        if current_id is not None:
            base_id = current_id.rsplit("_", 1)[0]
            if base_id == locus_id:
                candidates.append("".join(current_seq))

    with open(nucl_orfs_file) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        flush()

    best_match = None
    best_diff = float('inf')
    for nt_seq in candidates:
        trim = len(nt_seq) % 3
        nt_trimmed = nt_seq[:-trim] if trim else nt_seq
        seq_obj = Seq(nt_trimmed)
        translated_len = len(str(seq_obj.translate(to_stop=True)))
        diff = abs(translated_len - target_protein_len)
        if diff < best_diff and diff <= LENGTH_TOLERANCE_AA:
            best_diff = diff
            best_match = nt_trimmed

    if best_match is not None:
        exact_len = target_protein_len * 3
        if len(best_match) >= exact_len:
            return best_match[:exact_len]
        return best_match
    return None


def main():
    if len(sys.argv) != 2:
        print("Usage: python3 prepare_hyphy_input_en.py OGxxxxxxx")
        sys.exit(1)

    og_id = sys.argv[1]
    msa_file = find_msa_file(og_id)
    if not msa_file:
        print(f"ERROR: cannot find the alignment for {og_id}")
        sys.exit(1)

    print(f"Using alignment: {msa_file}")

    protein_records = list(SeqIO.parse(msa_file, "fasta"))
    print(f"Sequences in the group: {len(protein_records)}")

    out_dir = os.path.join(PROJECT_DIR, "05_selection_hyphy", og_id)
    os.makedirs(out_dir, exist_ok=True)

    codon_records = []
    n_ok, n_fail = 0, 0

    for rec in protein_records:
        parsed = parse_header(rec.id)
        if not parsed:
            print(f"  [!] Could not parse header: {rec.id}")
            n_fail += 1
            continue
        short_genome, locus_id, method = parsed

        genome_folder = find_genome_folder(short_genome)
        if not genome_folder:
            print(f"  [!] Cannot find folder for {short_genome}")
            n_fail += 1
            continue

        aligned_protein = str(rec.seq)
        ungapped_len = len(aligned_protein.replace("-", ""))

        nt_seq = None
        if method == "miniprot":
            nt_seq = get_nucleotide_from_miniprot(genome_folder, locus_id, ungapped_len)
        elif method == "getorf":
            nt_seq = get_nucleotide_from_getorf(genome_folder, locus_id, ungapped_len)

        if nt_seq is None:
            print(f"  [!] Could not recover nucleotides for {rec.id}")
            n_fail += 1
            continue

        try:
            from codon_thread_en import thread_codons
            codon_seq = thread_codons(aligned_protein, nt_seq)
        except ValueError as e:
            print(f"  [!] {rec.id}: {e}")
            n_fail += 1
            continue

        codon_records.append((rec.id, codon_seq))
        n_ok += 1

    out_fasta = os.path.join(out_dir, f"{og_id}_codon_alignment.fasta")
    with open(out_fasta, "w") as f:
        for name, seq in codon_records:
            f.write(f">{name}\n{seq}\n")

    print(f"\nDone: {n_ok} sequences with recovered codons, {n_fail} failed.")
    print(f"Codon alignment saved to: {out_fasta}")

    if n_fail > 0:
        print(f"\nNOTE: HyPhy needs the SAME set of taxa in the tree and")
        print(f"in the alignment. If you dropped sequences here, you may")
        print(f"need to prune the .treefile to remove those same tips")
        print(f"before running HyPhy.")


if __name__ == "__main__":
    main()
