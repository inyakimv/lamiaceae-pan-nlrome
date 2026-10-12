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
| EMBOSS (getorf) | 6.6.0.0 | `-find 1 -minsize 300` (protein ORFs); `-find 3 -minsize 300` (nucleotide recovery in `selection/`) |
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
| ete3 | 3.1.3 | species-tree plotting; tree pruning (`phylogeny/prune_tree_en.py`) |

## Pipeline and folders

| Folder | Purpose |
|---|---|
| `qc/` | Assembly statistics (SeqKit), gene-space completeness (BUSCO), Mash identity check and Table 2 |
| `annotation/` | NLR identification (NLR-Annotator), ANNEVO gene prediction, protein extraction, HMMER domain search and per-genome summary tables |
| `method_comparison/` | Comparison of NLR detection methods |
| `orthology/` | Orthogroup inference (OrthoFinder) and tandem-duplication analysis |
| `ploidy/` | Ploidy assignment per genome (CCDB, literature, assembly structure) |
| `figures/` | Scripts that generate the figures (white background, no gridlines) |
| `phylogeny/` | Scripts to build and prune the orthogroup maximum-likelihood trees (IQ-TREE; protein alignments from OrthoFinder) |
| `selection/` | Codon alignments and selection tests (HyPhy FEL/MEME) |
| `structure/` | Structural comparison of NLR candidates (AlphaFold, Foldseek) |
| `validation/` | Validation against the RefPlantNLR reference set (BLAST) |
| `data/` | Input and intermediate tables used by the scripts (CSV) |
| `results_tables/` | Result tables (e.g. Mash distances) |
| `output/figures/` | Source data behind each figure panel (CSV) |

## Scripts by results section

### 4.1 Size and tandem organization of the NLR repertoire across 46 Lamiaceae genomes

| Script | What it does | Output |
|---|---|---|
| `qc/qc_seqkit_busco_flat.sh` | SeqKit statistics and BUSCO for each assembly | `data/qc_summary.csv` |
| `qc/verify_seqkit_columns.sh` | Checks the SeqKit v2.13.0 column positions (N50, GC) | printed check |
| `qc/mash_mentha_check.sh` | Mash sketches and pairwise distances among *Mentha* assemblies and reference species | `results_tables/mash_mentha_check.tsv` |
| `qc/combine_qc_nlr_fixed.sh` | Merges QC statistics with NLR counts | merged QC and NLR table |
| `qc/qc_table2_and_correlations.py` | Assembly QC summary table, outlier fence and Pearson correlations (genome size vs NLR count and density) | printed tables |
| `ploidy/fetch_ploidy_ccdb.py` | Chromosome-number records from the Chromosome Counts Database | `data/ccdb_summary.csv` |
| `ploidy/assembly_chromosome_counts.py` | Sequence counts per assembly, including sequences of at least 1 Mb | `data/assembly_chromosome_counts.csv` |
| `ploidy/build_ploidy_worksheet.py` | Combines CCDB, assembly and BUSCO evidence into one worksheet | `data/ploidy_worksheet.csv` |
| `ploidy/make_binary_ploidy.py` | Diploid / polyploid table | `data/ploidy_binary.csv` |
| `figures/plot_nlr_ploidy.py` | NLR candidates per genome next to the species tree, coloured by ploidy | `output/figures/Fig1a_source_data.csv` |
| `orthology/run_orthofinder.sh` | Orthogroup inference with OrthoFinder (run `Results_Jul07`) | OrthoFinder results |
| `orthology/analyze_tandem_46species.py` | Tandem clusters of the most expanded NLR orthogroup per species (threshold 500 kb) | `data/tandem_analysis_46species.csv` |
| `orthology/tandem_summary_stats.py` | Median, range and correlations for the tandem analysis | printed summary |
| `figures/plot_tandem_hist.py` | Distribution of tandem-clustered loci across species | `output/figures/Fig1b_tandem_hist_source_data.csv` |

### 4.2 Cross-validation between NLR-Annotator and ANNEVO+HMMER

| Script | What it does | Output |
|---|---|---|
| `annotation/run_annevo_batch.sh` | Gene prediction with ANNEVO for each genome | `<genome>.annevo.gff3` |
| `annotation/extract_annevo_proteins.sh` | Protein sequences from the ANNEVO gene models (AGAT) | `<genome>.annevo_proteins.fasta` |
| `annotation/run_hmmer_nlr_search.sh` | `hmmscan` of the proteins against the NLR profiles (default parameters) | `<genome>.domtblout` |
| `annotation/nlr_profiles.hmm` | Seven Pfam profiles: NB-ARC, LRR_1, LRR_3, LRR_4, LRR_8, TIR, TIR_2 | profile database |
| `method_comparison/compare_nlr_methods_refined.py` | Per-species complete-like candidates (NB-ARC and at least one LRR in the same protein) next to NLR-Annotator counts | `data/nlr_method_comparison_detailed.csv` |
| `method_comparison/compare_nlr_loci.py` | Overlap of NLR-Annotator loci and ANNEVO+HMMER genes (same scaffold, at least one shared base, strand ignored) | `data/nlr_locus_comparison.csv` |
| `method_comparison/gene_size_distribution.py` | Locus length of each method, pooled over genomes; Mann-Whitney U test | `data/gene_lengths_nlr_annotator.csv`, `data/gene_lengths_annevo.csv` |
| `figures/plot_method_comparison.py` | Both panels of the comparison figure; Pearson and Spearman correlations, Cliff's delta | `output/figures/Fig2_method_comparison_source_data_a.csv`, `..._b_summary.csv` |

`method_comparison/compare_nlr_methods.py` is the first version of the comparison (it counted proteins with NB-ARC only) and is kept for reference.

### 4.3 Positive selection in five NLR orthogroups across ploidy levels and degrees of expansion concentration

| Script | What it does | Output |
|---|---|---|
| `selection/prepare_hyphy_input_en.py` | Recovers the in-frame nucleotide sequence of each protein (miniprot GFF or `getorf`) and threads it into the OrthoFinder protein alignment | `<OG>_codon_alignment.fasta` |
| `phylogeny/prune_tree_en.py` | Prunes the IQ-TREE gene tree to the sequences of the codon alignment (ete3) | `<OG>_pruned.treefile` |
| `selection/prune_tree_to_alignment.py` | Independent check: prunes the tree with Biopython and compares tip sets and patristic distances with the `_pruned.treefile` used for HyPhy (`--compare`) | `<OG>_pruned_check.treefile` |
| `selection/run_hyphy_en.sh` | HyPhy FEL and MEME on the codon alignment and pruned tree (default parameters) | `<OG>_FEL.json`, `<OG>_MEME.json` |
| `selection/summarize_hyphy_results_en.py` | Site counts and significant sites from the FEL and MEME output | `hyphy_summary.csv`, `hyphy_significant_sites.csv` |
| `selection/selection_summary.py` | Copies kept and lost, focal species, ploidy, informative codons and FEL/MEME counts per orthogroup | `data/selection_summary.csv`, `data/selection_summary_composition.csv` |
| `figures/plot_selection_summary.py` | Fig. 3: sites under selection per orthogroup | `output/figures/Fig3_selection_source_data.csv` |
| `figures/plot_og_tree.py` | Fig. 4: circular gene tree of OG0000000 coloured by species, with the per-species block table | `output/figures/Fig4_OG0000000_tree_species_blocks.csv`, `output/figures/Fig4_OG0000000_tree_tips.csv` |

The codon alignment, pruned tree and HyPhy output of each of the five orthogroups are in `data/selection/<OG>/`, and the tree of OG0000000 used for Fig. 4 is `data/OG0000000.treefile`.



### 4.4 Computational validation against RefPlantNLR

| Script | What it does | Output |
|---|---|---|
| `structure/extract_proteins_for_blast_en.py` | Ungapped protein sequences of the members of each orthogroup, from the `Results_Jul07` OrthoFinder alignments | `<OG>_proteins_for_blast.fasta` |
| `validation/run_blast_validation_en.sh` | Builds the RefPlantNLR database and runs `blastp` (E-value 1e-10, up to five targets per query) | `data/blast_results/<OG>_blast_results.tsv` |
| `validation/summarize_blast_results_en.py` | Best hit (highest bit score) per protein, reference groups, identity and coverage | `data/blast_validation_summary.csv`, `data/blast_validation_by_orthogroup.csv` |
| `validation/run_nbarc_phylogeny_en.sh` | Stages `prepare` (NB-ARC domains with `hmmscan`, MAFFT, column trimming), `tree` (IQ-TREE, `-m MFP -bb 1000`) and `placement` | `11_nbarc_tree/` (see `data/nbarc_tree/`) |
| `validation/nbarc_tools_en.py` | Helper called by the script above: joins queries and RefPlantNLR, extracts the NB-ARC domain, trims the alignment, and reports the smallest clade with reference proteins (and UFBoot) for every Lamiaceae protein | `nbarc_tree_placement.csv` |

The RefPlantNLR protein set (Kourelis et al., 2021; 415 sequences) is not redistributed here; its source is given in the thesis (Section 3.8, Annex A20).

The NB-ARC tree (Newick, with ultrafast-bootstrap values), the identifier map and the placement table are in `data/nbarc_tree/`. Alignments that contain RefPlantNLR sequences are not redistributed.

*(Tables for the next results sections are added here.)*

## Annex A21: structural comparison

- `scripts/structure/foldseek_statistics_unordered_en.py`: unordered pairs, Mann–Whitney U test, permutation test and leave-one-orthogroup-out analysis.
- `scripts/structure/plot_foldseek_violins_en.py`: Fig. 5.
- Outputs (tables and figures) are in `results/structure/`.
