#!/usr/bin/env python3
"""
Finds the best candidate for the 5th HyPhy group: a POLYPLOID
species with LOW concentration in its top orthogroup
(distributed expansion pattern, not concentrated tandem) --
the quadrant not yet covered in the 2x2 design.

Usage:
    cd ~/project
    python3 find_5th_candidate_en.py

Requires:
    - 03_orthofinder/proteomes/OrthoFinder/Results_*/Orthogroups/tandem_analysis_46species.csv
    - 09_combined_results/qc_nlr_combined.csv (or whatever you
      named your QC+NLR table with the ploidy_flag column)
"""
import pandas as pd
import glob
import os

PROJECT_DIR = os.path.expanduser("~/project")

tandem_files = glob.glob(os.path.join(
    PROJECT_DIR, "03_orthofinder/proteomes/OrthoFinder/Results_*/Orthogroups/tandem_analysis_46species.csv"
))
qc_files = glob.glob(os.path.join(PROJECT_DIR, "**/qc_nlr_combined*.csv"), recursive=True)

if not tandem_files:
    print("ERROR: cannot find tandem_analysis_46species.csv")
    exit(1)
if not qc_files:
    print("ERROR: cannot find a qc_nlr_combined*.csv file")
    print("Manually adjust the path in this script if it has a different name.")
    exit(1)

tandem = pd.read_csv(tandem_files[0])
qc = pd.read_csv(qc_files[0])

merged = tandem.merge(qc[["assembly_file", "busco_duplicated_pct"]],
                       left_on="species", right_on="assembly_file", how="left")

polyploids = merged[merged["busco_duplicated_pct"] > 60].sort_values("pct_concentration")

print("POLYPLOID species sorted by LOWEST concentration (best candidate first):")
print(polyploids[["species", "top_orthogroup", "pct_concentration",
                   "pct_genes_in_tandem_clusters", "busco_duplicated_pct"]].to_string(index=False))
