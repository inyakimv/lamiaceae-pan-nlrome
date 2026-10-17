# lamiaceae-pan-nlrome
Scripts for the pan-NLRome analysis of 46 Lamiaceae genomes (NLR identification, orthology, phylogeny, selection and structural analyses).
Conda environment specifications are provided in `env_specs/`.
## Software versions

| Tool | Version | Notes |
|---|---|---|
| GNU Wget | 1.21.2 | genome download from NCBI GenBank |
| SeqKit | 2.13.0 | `seqkit stats -a -T` |
| BUSCO | 5.7.1 | genome mode, lineage eudicots_odb10, default parameters |
| Mash | 2.3 | `-k 21 -s 1000000 -S 42` (Mentha identity check) |
| NLR-Annotator | 2.1b | default parameters |
| ANNEVO | 2.3.2 | `-l Magnoliopsida -t 24`, GPU; Python 3.10.20, PyTorch 2.1.0 (CUDA 12.1) |
| AGAT | 1.7.0 | `agat_sp_extract_sequences.pl -p` |
| HMMER (hmmscan) | 3.3.2 | default parameters, `--cpu 8`; NB-ARC domain search of Section 4.4 with `--cut_ga` |
| Pfam | 37.1 | NB-ARC, LRR_1/3/4/8, TIR, TIR_2 profiles |
| miniprot | 0.18-r281 | default parameters (`--gff`) |
| EMBOSS (getorf) | 6.6.0.0 | `-find 1 -minsize 300` (protein ORFs); `-find 3 -minsize 300` (nucleotide recovery in `scripts/selection/`) |
| Biopython | 1.87 | custom scripts |
| OrthoFinder | 3.1.5 | `orthofinder -f <proteomes> -t 32 -a 32` |
| FAMSA (OrthoFinder environment) | 2.4.1 | protein alignments run inside OrthoFinder (orthogroup alignments of Section 4.3) |
| DIAMOND (OrthoFinder environment) | 2.2.2 | sequence searches run inside OrthoFinder |
| FastTree (OrthoFinder environment) | 2.2.0 | gene trees run inside OrthoFinder (not used for Section 4.3) |
| IQ-TREE | 3.1.2 | `-m MFP -bb 1000 -nt 4` (gene trees of Section 4.3; NB-ARC tree of Section 4.4) |
| MAFFT | 7.525 | NB-ARC domain alignment of Section 4.4 (`--auto --thread 8`) |
| HyPhy | 2.5.100 (MP) | FEL / MEME, default parameters |
| BLAST+ (blastp, makeblastdb) | 2.12.0+ | `-evalue 1e-10 -max_target_seqs 5`, other parameters default (Section 4.4) |
| AlphaFold | 3 (AlphaFold Server) | |
| Foldseek | 10.941cd33 | |
| Python | 3.11.13 | custom scripts |
| pandas / NumPy / SciPy | 3.0.3 / 2.4.6 / 1.17.1 | tables, statistics |
| matplotlib | 3.11.1 | figures |
| ete3 | 3.1.3 | species-tree plotting; tree pruning (`scripts/phylogeny/prune_tree_en.py`) |
| MAFFT | 7.490 | clade-level alignments of Section 4.5 (`--auto --thread 8`) |
| IQ-TREE | 3.1.4 | clade-level trees of Section 4.5 (`-m MFP -mset LG,WAG,JTT,Q.plant,Q.pfam -mrate G4 -B 1000 -alrt 1000 -T 8`) |
| cd-hit | 4.8.1 | 85% identity clustering of NB-ARC domains (Section 4.5) |

## Pipeline and folders

| Folder | Purpose |
|---|---|
| `scripts/qc/` | Assembly statistics (SeqKit), gene-space completeness (BUSCO), Mash identity check and Table 2 |
| `scripts/annotation/` | NLR identification (NLR-Annotator), ANNEVO gene prediction, protein extraction, HMMER domain search and per-genome summary tables |
| `scripts/method_comparison/` | Comparison of NLR detection methods |
| `scripts/orthology/` | Orthogroup inference (OrthoFinder) and tandem-duplication analysis |
| `scripts/ploidy/` | Ploidy assignment per genome (CCDB, literature, assembly structure) |
| `scripts/figures/` | Scripts that generate the figures (white background, no gridlines) |
| `scripts/phylogeny/` | Scripts to build and prune the orthogroup maximum-likelihood trees (IQ-TREE; protein alignments from OrthoFinder) |
| `scripts/selection/` | Codon alignments and selection tests (HyPhy FEL/MEME) |
| `scripts/structure/` | Structural comparison of NLR candidates (AlphaFold, Foldseek) |
| `scripts/validation/` | Validation against the RefPlantNLR reference set (BLAST) |
| `data/` | Input tables (ploidy worksheet, CCDB, assembly counts), codon alignments and HyPhy output, BLAST results, NB-ARC tree, NLR HMM profiles (`data/hmm/`) |
| `results/` | Computed tables and figures by topic (`qc/`, `ploidy/`, `methods/`, `orthology/`, `structure/`, `figure_source_data/`) |
| `results/figure_source_data/` | Source data behind each figure panel (CSV) |
| `lamiaceae-nlr-conserved-core/` | Section 4.5: NB-ARC family tree, conserved NLR core of Lamiaceae (`scripts/conserved_core/`: scripts 01 to 07 and own README; `results/conserved_core/`: `figures/` and `tables/`) |

## Scripts by results section

### 4.1 Size and tandem organization of the NLR repertoire across 46 Lamiaceae genomes

| Script | What it does | Output |
|---|---|---|
| `scripts/qc/qc_seqkit_busco_flat.sh` | SeqKit statistics and BUSCO for each assembly | `results/qc/qc_summary.csv` |
| `scripts/qc/verify_seqkit_columns.sh` | Checks the SeqKit v2.13.0 column positions (N50, GC) | printed check |
| `scripts/qc/mash_mentha_check.sh` | Mash sketches and pairwise distances among *Mentha* assemblies and reference species | `results/qc/mash_mentha_check.tsv` |
| `scripts/qc/combine_qc_nlr_fixed.sh` | Merges QC statistics with NLR counts | merged QC and NLR table |
| `scripts/qc/qc_table2_and_correlations.py` | Assembly QC summary table, outlier fence and Pearson correlations (genome size vs NLR count and density) | printed tables |
| `scripts/ploidy/fetch_ploidy_ccdb.py` | Chromosome-number records from the Chromosome Counts Database | `data/ccdb_summary.csv` |
| `scripts/ploidy/assembly_chromosome_counts.py` | Sequence counts per assembly, including sequences of at least 1 Mb | `data/assembly_chromosome_counts.csv` |
| `scripts/ploidy/build_ploidy_worksheet.py` | Combines CCDB, assembly and BUSCO evidence into one worksheet | `data/ploidy_worksheet.csv` |
| `scripts/ploidy/make_binary_ploidy.py` | Diploid / polyploid table | `results/ploidy/ploidy_binary.csv` |
| `scripts/figures/plot_nlr_ploidy.py` | NLR candidates per genome next to the species tree, coloured by ploidy | `results/figure_source_data/Fig1a_source_data.csv` |
| `scripts/orthology/run_orthofinder.sh` | Orthogroup inference with OrthoFinder (run `Results_Jul07`) | OrthoFinder results |
| `scripts/orthology/analyze_tandem_46species.py` | Tandem clusters of the most expanded NLR orthogroup per species (threshold 500 kb) | `results/orthology/tandem_analysis_46species.csv` |
| `scripts/orthology/tandem_summary_stats.py` | Median, range and correlations for the tandem analysis | printed summary |
| `scripts/figures/plot_tandem_hist.py` | Distribution of tandem-clustered loci across species | `results/figure_source_data/Fig1b_tandem_hist_source_data.csv` |

### 4.2 Cross-validation between NLR-Annotator and ANNEVO+HMMER

| Script | What it does | Output |
|---|---|---|
| `scripts/annotation/run_annevo_batch.sh` | Gene prediction with ANNEVO for each genome | `<genome>.annevo.gff3` |
| `scripts/annotation/extract_annevo_proteins.sh` | Protein sequences from the ANNEVO gene models (AGAT) | `<genome>.annevo_proteins.fasta` |
| `scripts/annotation/run_hmmer_nlr_search.sh` | `hmmscan` of the proteins against the NLR profiles (default parameters) | `<genome>.domtblout` |
| `data/hmm/nlr_profiles.hmm` | Seven Pfam profiles: NB-ARC, LRR_1, LRR_3, LRR_4, LRR_8, TIR, TIR_2 | profile database |
| `scripts/method_comparison/compare_nlr_methods_refined.py` | Per-species complete-like candidates (NB-ARC and at least one LRR in the same protein) next to NLR-Annotator counts | `results/methods/nlr_method_comparison_detailed.csv` |
| `scripts/method_comparison/compare_nlr_loci.py` | Overlap of NLR-Annotator loci and ANNEVO+HMMER genes (same scaffold, at least one shared base, strand ignored) | `results/methods/nlr_locus_comparison.csv` |
| `scripts/method_comparison/gene_size_distribution.py` | Locus length of each method, pooled over genomes; Mann-Whitney U test | `results/methods/gene_lengths_nlr_annotator.csv`, `results/methods/gene_lengths_annevo.csv` |
| `scripts/figures/plot_method_comparison.py` | Both panels of the comparison figure; Pearson and Spearman correlations, Cliff's delta | `results/figure_source_data/Fig2_method_comparison_source_data_a.csv`, `..._b_summary.csv` |

`scripts/method_comparison/compare_nlr_methods.py` is the first version of the comparison (it counted proteins with NB-ARC only) and is kept for reference.

### 4.3 Positive selection in five NLR orthogroups across ploidy levels and degrees of expansion concentration

| Script | What it does | Output |
|---|---|---|
| `scripts/selection/prepare_hyphy_input_en.py` | Recovers the in-frame nucleotide sequence of each protein (miniprot GFF or `getorf`) and threads it into the OrthoFinder protein alignment | `<OG>_codon_alignment.fasta` |
| `scripts/phylogeny/prune_tree_en.py` | Prunes the IQ-TREE gene tree to the sequences of the codon alignment (ete3) | `<OG>_pruned.treefile` |
| `scripts/selection/prune_tree_to_alignment.py` | Independent check: prunes the tree with Biopython and compares tip sets and patristic distances with the `_pruned.treefile` used for HyPhy (`--compare`) | `<OG>_pruned_check.treefile` |
| `scripts/selection/run_hyphy_en.sh` | HyPhy FEL and MEME on the codon alignment and pruned tree (default parameters) | `<OG>_FEL.json`, `<OG>_MEME.json` |
| `scripts/selection/summarize_hyphy_results_en.py` | Site counts and significant sites from the FEL and MEME output | `hyphy_summary.csv`, `hyphy_significant_sites.csv` |
| `scripts/selection/selection_summary.py` | Copies kept and lost, focal species, ploidy, informative codons and FEL/MEME counts per orthogroup | `data/selection_summary.csv`, `data/selection_summary_composition.csv` |
| `scripts/figures/plot_selection_summary.py` | Fig. 3: sites under selection per orthogroup | `results/figure_source_data/Fig3_selection_source_data.csv` |
| `scripts/figures/plot_og_tree.py` | Fig. 4: circular gene tree of OG0000000 coloured by species, with the per-species block table | `results/figure_source_data/Fig4_OG0000000_tree_species_blocks.csv`, `results/figure_source_data/Fig4_OG0000000_tree_tips.csv` |

The codon alignment, pruned tree and HyPhy output of each of the five orthogroups are in `data/selection/<OG>/`, and the tree of OG0000000 used for Fig. 4 is `data/OG0000000.treefile`.

### 4.4 Computational validation against RefPlantNLR

| Script | What it does | Output |
|---|---|---|
| `scripts/structure/extract_proteins_for_blast_en.py` | Ungapped protein sequences of the members of each orthogroup, from the `Results_Jul07` OrthoFinder alignments | `<OG>_proteins_for_blast.fasta` |
| `scripts/validation/run_blast_validation_en.sh` | Builds the RefPlantNLR database and runs `blastp` (E-value 1e-10, up to five targets per query) | `data/blast_results/<OG>_blast_results.tsv` |
| `scripts/validation/summarize_blast_results_en.py` | Best hit (highest bit score) per protein, reference groups, identity and coverage | `data/blast_validation_summary.csv`, `data/blast_validation_by_orthogroup.csv` |
| `scripts/validation/run_nbarc_phylogeny_en.sh` | Stages `prepare` (NB-ARC domains with `hmmscan`, MAFFT, column trimming), `tree` (IQ-TREE, `-m MFP -bb 1000`) and `placement` | `11_nbarc_tree/` (see `data/nbarc_tree/`) |
| `scripts/validation/nbarc_tools_en.py` | Helper called by the script above: joins queries and RefPlantNLR, extracts the NB-ARC domain, trims the alignment, and reports the smallest clade with reference proteins (and UFBoot) for every Lamiaceae protein | `nbarc_tree_placement.csv` |

The RefPlantNLR protein set (Kourelis et al., 2021; 415 sequences) is not redistributed here; its source is given in the thesis (Section 3.8, Annex A20).

The NB-ARC tree (Newick, with ultrafast-bootstrap values), the identifier map and the placement table are in `data/nbarc_tree/`. Alignments that contain RefPlantNLR sequences are not redistributed.

### 4.5 Conserved NLR core of Lamiaceae

| Script | What it does | Output |
|---|---|---|
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/01_partition_refs.py` | Cuts the NB-ARC family tree (Lamiaceae, outgroups and RefPlantNLR anchors) into clades of at most 300 leaves | `refs_collapsed.treefile`, `refs_clade_info.csv`, `refs_clade_members.tsv` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/02_clade_stats.py` | Expands cd-hit (85%) clusters to genes; species present and copies per species for each clade | `stats_clade_conservation.tsv`, `stats_clade_genes.tsv` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/03_fig_core_clades_circular.py` | Circular figure of the 85 clades, rooted on the TNL/RNL clade; shades conserved clades (at least 35 species, at most 3 genes per species; darker if the median is one copy) | `lamiaceae-nlr-conserved-core/results/conserved_core/figures/Fig_nbarc_core_clades_circular_rooted.png/.pdf` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/04_inspect_clade.py` | Members, copies per species, tandem arrays and neighbouring clades of one clade | `<clade>_genes.tsv`, `<clade>_nbarc.fasta` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/05_mada_check.py` | N-terminal (MADA-like) check of the proteins of one clade against ZAR1 and NRC references | `<clade>_mada.tsv` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/06_clade_fasta.py` | FASTA of a clade plus selected references and outgroups | `<clade>_tree.fasta` |
| `lamiaceae-nlr-conserved-core/scripts/conserved_core/07_fig_clade_tree.py` | Rectangular clade tree (collapsed or full) from the IQ-TREE output | `lamiaceae-nlr-conserved-core/results/conserved_core/figures/Fig_S_Clade18_tree_*.png/.pdf` |

Tables in `lamiaceae-nlr-conserved-core/results/conserved_core/tables/`: per-species NLR summary (NB-ARC counts, ADR1, NRG1, RPW8, TIR, EDS1, PAD4, SAG101, tblastn hits, copies in the four conserved clades and ZAR1 tblastn hits), sensitivity of the conserved-clade definition to the thresholds, and host species of the reference NLRs by clade.

## Annex A21: structural comparison

- `scripts/structure/foldseek_statistics_unordered_en.py`: unordered pairs, Mann–Whitney U test, permutation test and leave-one-orthogroup-out analysis.
- `scripts/structure/plot_foldseek_violins_en.py`: Fig. 5.
- Outputs (tables and figures) are in `results/structure/`.

## Notes on the scripts
- Paths (/data/tfg_users/mollai/...) are those of the server where the analysis was run. To reproduce, change the variables 
- `scripts/archive_not_used_in_thesis/`: earlier or one-off versions kept for the record.
