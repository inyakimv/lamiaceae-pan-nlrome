# lamiaceae-pan-nlrome
Scripts for the pan-NLRome analysis of 46 Lamiaceae genomes (NLR identification, orthology, phylogeny, selection and structural analyses).
## Software versions

| Tool | Version | Notes |
|---|---|---|
| GNU Wget | 1.21.2 | genome download from NCBI GenBank |
| SeqKit | 2.13.0 | `seqkit stats -a -T` |
| BUSCO | 5.7.1 | genome mode, lineage eudicots_odb10, default parameters |
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
