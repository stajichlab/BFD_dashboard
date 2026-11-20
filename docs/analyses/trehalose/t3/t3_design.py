#!/usr/bin/env python3
"""T3 step 1: design check. How many genera have both labelled and unlabelled kept species? Touches NO outcome metric.

Usage: t3_design.py <J13 data dir (th_genomes.csv, th_sp.csv, th_univcov.csv, th_univ.csv)> <out dir>
Writes design_genera.tsv (per genus counts) and labelled_species.tsv (every labelled species with the epithet that matched).
"""
import re
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
import t3_lib as L

data, out = Path(sys.argv[1]), Path(sys.argv[2])
out.mkdir(parents=True, exist_ok=True)
X = L.load_species(data)
T = L.design_table(X)
T.to_csv(out / "design_genera.tsv", sep="\t", index=False)
lab = X[X.label][["GENUS", "SPECIES", "CLASS", "keep"]].copy()
lab["matched"] = lab.SPECIES.map(lambda s: (re.search(L.DRY_EPITHET, s) or [None])[0])
lab.sort_values(["GENUS", "SPECIES"]).to_csv(out / "labelled_species.tsv", sep="\t", index=False)
print(f"species representatives: {len(X)}; kept: {int(X.keep.sum())}; labelled (any quality): {int(X.label.sum())}; labelled and kept: {int((X.label & X.keep).sum())}")
print(f"genera with >=1 labelled kept species: {int((T.n_label >= 1).sum())}; eligible (>= {L.MIN_FLAG} labelled and >= {L.MIN_NON} unlabelled kept): {int(T.eligible.sum())}")
print(T[T.n_label >= 1].to_string(index=False))
