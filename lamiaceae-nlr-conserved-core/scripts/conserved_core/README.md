# NB-ARC family tree of Lamiaceae NLRs: clade statistics, conserved core and figures
Scripts to partition the NB-ARC phylogeny of Lamiaceae NLRs into clades, annotate them with RefPlantNLR reference NLRs, compute per-clade conservation across 46 Lamiaceae genomes, and draw the circular overview figure used in the paper.

## Overview

| Step | Script | Input | Output |
|---|---|---|---|
| 1 | `01_partition_refs.py` | NB-ARC tree (Lamiaceae + *Arabidopsis thaliana* + *Klebsormidium nitens* + RefPlantNLR anchors), RefPlantNLR Table S1 | `refs_collapsed.treefile`, `refs_clade_info.csv`, `refs_clade_members.tsv`, `refs_anchors.tsv` |
| 2 | `02_clade_stats.py` | cd-hit (85%) `.clstr`, step 1 members and info | `<prefix>_clade_conservation.tsv`, `<prefix>_clade_genes.tsv` |
| 3 | `03_fig_core_clades_circular.py` | step 1 folder, step 2 tables | `<out>.png`, `<out>.pdf` |
| 4 | `04_inspect_clade.py` | one clade id, step 1 folder, cd-hit `.clstr` | `<prefix>_genes.tsv` (members, copies per species, contig and gene index), `<prefix>_nbarc.fasta`; prints tandem arrays and neighbouring clades |
| 5 | `05_mada_check.py` | step 4 gene table, protein FASTAs, optional RefPlantNLR FASTA | N-terminal (MADA-like) table; prints the N-termini of ZAR1 and NRC references |
| 6 | `06_clade_fasta.py` | step 4 gene table, protein FASTAs, RefPlantNLR FASTA | FASTA of a clade plus selected references and outgroups, ready for MAFFT and IQ-TREE |
| 7 | `07_fig_clade_tree.py` | IQ-TREE `.treefile` of a clade | rectangular tree figure (`collapsed` or `full`) |

## Method in brief

1. **Clades.** The tree is midpoint-rooted only to cut it into maximal subtrees of at most 300 leaves (85 clades for the published tree); the figure is then re-rooted (step 3). Each clade is annotated with the RefPlantNLR references it contains and their class (CNL, TNL, RNL, from Table S1).
2. **Gene counts.** The tree was built from cd-hit representatives (85% identity). Each representative is expanded with the members of its cluster, so clade sizes in steps 2 and 3 are counted in genes, not tree leaves.
3. **Conserved low-copy clades.** A clade is flagged when it contains Lamiaceae genes from at least 35 of the 46 species and has at most 3 genes per species on average (thresholds in `core()` in script 03; the conserved clades were stable across 30–45 species and ≤2–5 copies per species, see `Table_S_threshold_sensitivity.tsv`). Clades whose median is one copy per species are drawn darker (`stable()`), the rest are widespread but locally expanded.
4. **Figure.** Circular layout, inside to outside: collapsed tree; share of leaves by origin (Lamiaceae, Arabidopsis, reference CNL, TNL, RNL); number of Lamiaceae species present (grey scale); Lamiaceae genes per clade (log scale).

## Usage (run from scripts/conserved_core/)

```bash
python3 01_partition_refs.py nbarc_hmm_plus_refs.treefile TableS1_description.xlsx out/refs 300
python3 02_clade_stats.py nbarc_domains_cdhit85.fasta.clstr out/refs_clade_members.tsv out/refs_clade_info.csv out/stats
python3 03_fig_core_clades_circular.py out out/stats_clade_conservation.tsv out/stats_clade_genes.tsv out/Fig_nbarc_core_clades_circular Clade18
```

Clade-level analyses (steps 4 to 7), shown for Clade 18 (ADR1/NRG1):

```bash
python3 04_inspect_clade.py Clade18 out nbarc_domains_cdhit85.fasta.clstr out/Clade18 nbarc_domains_combined.fasta
python3 06_clade_fasta.py out/Clade18_genes.tsv annevo_proteins RefPlantNLR_proteins.fasta 'ADR1|NRG1' out/Clade18_tree.fasta '^ZAR1$|^NRC0$'
mafft --auto --thread 8 out/Clade18_tree.fasta > out/Clade18_aln.fasta
iqtree3 -s out/Clade18_aln.fasta -m MFP -mset LG,WAG,JTT,Q.plant,Q.pfam -mrate G4 -B 1000 -alrt 1000 -T 8 -o OUT_ZAR1,OUT_NRC0 --prefix out/Clade18_tree
python3 07_fig_clade_tree.py out/Clade18_tree.treefile out/Fig_S_Clade18_tree_collapsed collapsed
```

N-terminal check of a clade (for example Clade 76 against NRC references):

```bash
python3 05_mada_check.py out/Clade76_genes.tsv annevo_proteins out/Clade76_mada.tsv RefPlantNLR_proteins.fasta
```

The optional last argument of step 3 sets the root: `midpoint` (display only) or a clade used as outgroup. `Clade18` (the TNL/RNL clade, which also contains the *Klebsormidium* sequence) roots the tree between TNL/RNL and the CNL clades, which is the root supported by the full NB-ARC phylogeny.

The first argument of script 03 is the folder holding the step 1 outputs (prefix `refs`).

Paths in the examples (`out/`, `annevo_proteins/`, `RefPlantNLR_proteins.fasta`) are those of the analysis server and must be adapted.

## Requirements

Python 3 with `pandas`, `matplotlib`, `biopython` and `openpyxl`. Steps 4 to 7 additionally use MAFFT and IQ-TREE 3 (run in a separate conda environment if matplotlib is not installed there; the figure scripts need matplotlib).

## Input conventions

- Leaf names follow `<genome>_genomic|<contig>-g<N>.t1` (ANNEVO gene models); references are prefixed `REF|`.
- Cluster members in the `.clstr` file use the same names, so representatives can be matched to clades.

## Citing

If you use these scripts, please cite the associated paper (reference to be added) and the software they depend on: matplotlib, Biopython, pandas, cd-hit, and the RefPlantNLR resource.

## AI assistance

Parts of this code were drafted with an AI assistant (Claude, Anthropic) and were run, reviewed and validated by the authors. Edit or remove this section according to the policy of the target journal.
