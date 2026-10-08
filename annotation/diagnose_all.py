#!/usr/bin/env python3
"""
Script de diagnostico: para OG0000038, muestra el 'diff' real
(diferencia en aminoacidos) del MEJOR candidato encontrado para
CADA secuencia, incluso las que fallan el filtro de <=10 aa.
Esto permite clasificar los fallos en "pequenos" (arreglables
ampliando el margen) vs "grandes" (requieren revision manual).
"""
import sys, os, glob, re, subprocess
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq

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
    parts = header.split("|")
    if len(parts) != 3:
        return None
    short_genome_field, locus_id, method = parts
    m = re.search(r'(GC[AF]_\d+\.\d+)', short_genome_field)
    if not m:
        return None
    return m.group(1), locus_id, method

def find_genome_folder(short_genome):
    if not os.path.isdir(REFINED_DIR):
        return None
    for name in os.listdir(REFINED_DIR):
        if name.startswith(short_genome):
            return os.path.join(REFINED_DIR, name)
    return None

def best_diff_miniprot(genome_folder, locus_id, target_protein_len):
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
                mrna_strand[attr_dict.get("ID", "")] = (seqid, strand)
            elif feature == "CDS":
                attr_dict = dict(a.split("=", 1) for a in attrs.split(";") if "=" in a)
                parent = attr_dict.get("Parent", "")
                ph = int(phase) if phase.isdigit() else 0
                cds_by_mrna[parent].append((int(start), int(end), ph))
    best = None
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
        translated_len = len(str(seq_obj.translate(to_stop=True)))
        diff = abs(translated_len - target_protein_len)
        if best is None or diff < best:
            best = diff
    return best

def best_diff_getorf(genome_folder, locus_id, target_protein_len):
    nucl_orfs_file = os.path.join(genome_folder, "nlr_plain_orfs_nucl.fasta")
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
    best = None
    for nt_seq in candidates:
        trim = len(nt_seq) % 3
        seq_obj = Seq(nt_seq[:-trim] if trim else nt_seq)
        translated_len = len(str(seq_obj.translate(to_stop=True)))
        diff = abs(translated_len - target_protein_len)
        if best is None or diff < best:
            best = diff
    return best

def main():
    og_id = sys.argv[1] if len(sys.argv) > 1 else "OG0000038"
    msa_file = find_msa_file(og_id)
    records = list(SeqIO.parse(msa_file, "fasta"))

    results = []
    for rec in records:
        parsed = parse_header(rec.id)
        if not parsed:
            results.append((rec.id, "parse_error", None))
            continue
        short_genome, locus_id, method = parsed
        genome_folder = find_genome_folder(short_genome)
        if not genome_folder:
            results.append((rec.id, "no_folder", None))
            continue
        ungapped_len = len(str(rec.seq).replace("-", ""))
        if method == "miniprot":
            diff = best_diff_miniprot(genome_folder, locus_id, ungapped_len)
        else:
            diff = best_diff_getorf(genome_folder, locus_id, ungapped_len)
        results.append((rec.id.split("|")[1] + "|" + method, "ok", diff))

    # Ordenar por diff descendente (peores primero), None al final
    results.sort(key=lambda x: (x[2] is None, -(x[2] or 0)))

    print(f"{'locus|method':<45} {'status':<10} {'diff (aa)'}")
    for name, status, diff in results:
        print(f"{name:<45} {status:<10} {diff}")

    diffs = [d for _, s, d in results if d is not None]
    print(f"\nTotal: {len(results)}")
    print(f"Con diff <= 10 (ya se recuperan): {sum(1 for d in diffs if d <= 10)}")
    print(f"Con diff 11-20 (posiblemente arreglable): {sum(1 for d in diffs if 10 < d <= 20)}")
    print(f"Con diff > 20 (probable estructura de intrones real): {sum(1 for d in diffs if d > 20)}")
    print(f"Sin ningun candidato encontrado: {sum(1 for _,s,d in results if d is None)}")

if __name__ == "__main__":
    main()
