#!/usr/bin/env python3
"""
Collapse ploidy levels (2x, 4x, ...) into "diploid" / "polyploid" for plot_nlr_ploidy.py.
Usage: python3 make_binary_ploidy.py ploidy.csv ploidy_binary.csv
"""
import sys

import pandas as pd

src, dst = sys.argv[1], sys.argv[2]
p = pd.read_csv(src)
p["ploidy"] = p["ploidy"].map(lambda v: "diploid" if str(v).strip() == "2x" else "polyploid")
p.to_csv(dst, index=False)
print(p["ploidy"].value_counts().to_string())
