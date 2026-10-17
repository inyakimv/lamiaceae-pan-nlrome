## Additions to the root README.md (flat folder layout)

### 1. Software versions (add these rows to the table)

| Tool | Version | Notes |
|---|---|---|
| MAFFT | 7.490 | clade-level alignments of Section 4.5 (`--auto --thread 8`) |
| IQ-TREE | 3.1.4 | clade-level trees of Section 4.5 (`-m MFP -mset LG,WAG,JTT,Q.plant,Q.pfam -mrate G4 -B 1000 -alrt 1000 -T 8`) |
| cd-hit | 4.8.1 | 85% identity clustering of NB-ARC domains (Section 4.5) |

### 2. Folders (add this row to the table)

| Folder | Purpose |
|---|---|
| `lamiaceae-nlr-conserved-core/` | Section 4.5: NB-ARC family tree, conserved NLR core of Lamiaceae (scripts 01 to 07, `figures/`, `tables/`, own README) |

### 3. New section (paste after Section 4.4)

#### 4.5 Conserved NLR core of Lamiaceae
| Script | What it does | Output |
|---|---|---|
| `lamiaceae-nlr-conserved-core/01_partition_refs.py` | Cuts the NB-ARC family tree (Lamiaceae, outgroups and RefPlantNLR anchors) into clades of at most 300 leaves | `refs_collapsed.treefile`, `refs_clade_info.csv`, `refs_clade_members.tsv` |
| `lamiaceae-nlr-conserved-core/02_clade_stats.py` | Expands cd-hit (85%) clusters to genes; species present and copies per species for each clade | `stats_clade_conservation.tsv`, `stats_clade_genes.tsv` |
| `lamiaceae-nlr-conserved-core/03_fig_core_clades_circular.py` | Circular figure of the 85 clades, rooted on the TNL/RNL clade; shades conserved clades (at least 35 species, at most 3 genes per species; darker if the median is one copy) | `lamiaceae-nlr-conserved-core/figures/Fig_nbarc_core_clades_circular_rooted.png/.pdf` |
| `lamiaceae-nlr-conserved-core/04_inspect_clade.py` | Members, copies per species, tandem arrays and neighbouring clades of one clade | `<clade>_genes.tsv`, `<clade>_nbarc.fasta` |
| `lamiaceae-nlr-conserved-core/05_mada_check.py` | N-terminal (MADA-like) check of the proteins of one clade against ZAR1 and NRC references | `<clade>_mada.tsv` |
| `lamiaceae-nlr-conserved-core/06_clade_fasta.py` | FASTA of a clade plus selected references and outgroups | `<clade>_tree.fasta` |
| `lamiaceae-nlr-conserved-core/07_fig_clade_tree.py` | Rectangular clade tree (collapsed or full) from the IQ-TREE output | `lamiaceae-nlr-conserved-core/figures/Fig_S_Clade18_tree_*.png/.pdf` |

Tables in `lamiaceae-nlr-conserved-core/tables/`: per-species NLR summary (NB-ARC counts, ADR1, NRG1, RPW8, TIR, EDS1, PAD4, SAG101, tblastn hits, copies in the four conserved clades and ZAR1 tblastn hits), sensitivity of the conserved-clade definition to the thresholds, and host species of the reference NLRs by clade.
