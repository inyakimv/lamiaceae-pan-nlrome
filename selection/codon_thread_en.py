"""
Core codon threading logic (similar to PAL2NAL):
takes an ALREADY ALIGNED protein sequence (with '-' gaps) and
the corresponding UNALIGNED nucleotide sequence (in frame,
length = 3 * ungapped protein length), and produces the
ALIGNED nucleotide version, inserting '---' wherever the
protein had a gap.
"""

def thread_codons(aligned_protein, nucleotide_seq):
    """
    aligned_protein: string, e.g. "MK-LST" (with '-' gaps)
    nucleotide_seq: string of gap-free nucleotides, length a
                    multiple of 3 (length = 3 * number of amino
                    acids NOT counting gaps), already in the
                    correct reading frame
    Returns: gapped nucleotide string ('---' for each amino
             acid gap), same aligned length as the protein.
    """
    ungapped_len = len(aligned_protein.replace("-", ""))
    expected_nt_len = ungapped_len * 3
    if len(nucleotide_seq) != expected_nt_len:
        raise ValueError(
            f"Length mismatch: ungapped protein={ungapped_len} aa "
            f"(expected {expected_nt_len} nt), actual nucleotides={len(nucleotide_seq)} nt"
        )

    result = []
    nt_pos = 0
    for aa in aligned_protein:
        if aa == "-":
            result.append("---")
        else:
            result.append(nucleotide_seq[nt_pos:nt_pos+3])
            nt_pos += 3
    return "".join(result)
