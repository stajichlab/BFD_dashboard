#!/usr/bin/env python3
"""T2: join FungalTraits 1.2 lifestyle labels to the candidate-loss genomes (genus level) and apply the pre-registered rule.

Usage: 05_lifestyle.py CANDIDATE_LOSSES_TSV FUNGALTRAITS_CSV OUT_PREFIX
Rule (README): T2 is killed if a candidate-loss genome has a FungalTraits primary lifestyle that is a free-living saprotroph
(any primary_lifestyle ending in '_saprotroph'). Labels are genus-level and are copied verbatim with their source; a genus that
FungalTraits lacks is reported as 'no label', never guessed. The FungalTraits primary_lifestyle column is read as the
`primary_lifestyle` field.
"""
import sys

import pandas as pd

SAPROTROPH_SUFFIX = "_saprotroph"


def label_genomes(cand, ft):
    """cand: DataFrame with column genus; ft: FungalTraits DataFrame with GENUS and primary_lifestyle (and Secondary_lifestyle)."""
    key = ft.rename(columns=str.strip).copy()
    gcol = next(c for c in key.columns if c.lower() in ("genus", "genus_"))
    key = key.drop_duplicates(gcol).set_index(gcol)
    out = cand.copy()
    for col in ("primary_lifestyle", "Secondary_lifestyle"):
        out[col] = out["genus"].map(key[col]) if col in key.columns else None
    out["primary_lifestyle"] = out["primary_lifestyle"].fillna("no label")
    out["free_living_saprotroph"] = out["primary_lifestyle"].str.endswith(SAPROTROPH_SUFFIX)
    return out


def t2_verdict(lab):
    n_sap = int(lab["free_living_saprotroph"].sum())
    if n_sap:
        genera = ", ".join(sorted(lab.loc[lab["free_living_saprotroph"], "genus"].unique()))
        return (f"T2 KILLED by its pre-registered rule: {n_sap} of {len(lab)} candidate-loss genomes are in genera whose "
                f"FungalTraits primary lifestyle is a free-living saprotroph ({genera})")
    return f"T2 not killed: no candidate-loss genome has a free-living saprotroph primary lifestyle ({len(lab)} genomes)"


def main(cand_tsv, ft_csv, prefix):
    cand = pd.read_csv(cand_tsv, sep="\t")
    ft = pd.read_csv(ft_csv, encoding="utf-8")
    lab = label_genomes(cand, ft)
    lab.to_csv(f"{prefix}_candidate_loss_lifestyle.tsv", sep="\t", index=False)
    by = lab.groupby(["genus", "primary_lifestyle", "free_living_saprotroph"]).size().reset_index(name="genomes")
    by.to_csv(f"{prefix}_lifestyle_by_genus.tsv", sep="\t", index=False)
    text = t2_verdict(lab) + "\n\n" + by.to_string(index=False) + "\n"
    open(f"{prefix}_t2_verdict.txt", "w").write(text)
    print(text)


if __name__ == "__main__":
    main(*sys.argv[1:4])
