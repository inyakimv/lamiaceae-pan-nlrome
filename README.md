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
| miniprot | 0.18-r281 | default parameters (`--gff`) |
| EMBOSS (getorf) | 6.6.0.0 | `-find 1 -minsize 300` |
| Biopython | 1.87 | custom scripts |
| OrthoFinder | 3.1.5 | `orthofinder -f <proteomes> -t 32 -a 32` |
| MAFFT | 7.490 (2021/Oct/30) | |
| IQ-TREE | 3.1.2 | |
| HyPhy | 2.5.100 (MP) | FEL / MEME |
| AlphaFold | 3 (AlphaFold Server) | |
| Foldseek | 10.941cd33 | |
| Python | 3.11.13 | custom scripts |
| pandas / NumPy / SciPy | 3.0.3 / 2.4.6 / 1.17.1 | tables, statistics |
| matplotlib | 3.11.1 | figures |
| ete3 | 3.1.3 | species-tree plotting |

## Pipeline and folders

| Folder | Purpose |
|---|---|
| `qc/` | Assembly statistics (SeqKit), gene-space completeness (BUSCO), Mash identity check and Table 2 |
| `annotation/` | NLR identification (NLR-Annotator), protein extraction (miniprot, getorf), ANNEVO-based annotation, per-genome summary tables |
| `method_comparison/` | Comparison of NLR detection methods |
| `orthology/` | Orthogroup inference (OrthoFinder) and tandem-duplication analysis |
| `ploidy/` | Ploidy assignment per genome (CCDB, literature, assembly structure) |
| `figures/` | Scripts that generate the figures (white background, no gridlines) |
| `phylogeny/` | Alignments and maximum-likelihood trees (MAFFT, IQ-TREE) |
| `selection/` | Codon alignments and selection tests (HyPhy FEL/MEME) |
| `structure/` | Structural comparison of NLR candidates (AlphaFold, Foldseek) |
| `validation/` | Validation against the RefPlantNLR reference set (BLAST) |
| `data/` | Input and intermediate tables used by the scripts (CSV) |
| `results_tables/` | Result tables (e.g. Mash distances) |
| `output/figures/` | Source data behind each figure panel (CSV) |

## Scripts for Section 4.1 and their Annex entries

| Script | What it does | Annex |
|---|---|---|
| `qc/qc_seqkit_busco_flat.sh` | SeqKit statistics and BUSCO for each assembly | A13 |
| `qc/verify_seqkit_columns.sh` | Checks the SeqKit v2.13.0 column positions (N50, GC) | A13 |
| `qc/mash_mentha_check.sh` | Mash sketches and pairwise distances among *Mentha* assemblies | A14 |
| `qc/qc_table2_and_correlations.py` | Table 2, outlier fence and Pearson correlations | A15 |
| `qc/combine_qc_nlr_fixed.sh` | Merges QC statistics with NLR counts (`data/qc_summary.csv`) | A15 |
| `ploidy/fetch_ploidy_ccdb.py` | Chromosome-number records from CCDB | A16 |
| `ploidy/assembly_chromosome_counts.py` | Sequence counts per assembly (including sequences of at least 1 Mb) | A16 |
| `ploidy/build_ploidy_worksheet.py` | Combines the evidence into one worksheet | A16 |
| `ploidy/make_binary_ploidy.py` | Diploid / polyploid table (`data/ploidy_binary.csv`) | A16 |
| `figures/plot_nlr_ploidy.py` | Fig. 1a, NLR candidates per genome on the species tree, coloured by ploidy | A16 |
| `orthology/analyze_tandem_46species.py` | Tandem clusters of the most expanded NLR orthogroup per species (threshold 500 kb) | A17 |
| `orthology/tandem_summary_stats.py` | Median, range and correlations for the tandem analysis | A17 |
| `figures/plot_tandem_hist.py` | Fig. 1b, distribution of tandem-clustered loci | A17 |

`orthology/run_orthofinder.sh` runs OrthoFinder (run `Results_Jul07`).
