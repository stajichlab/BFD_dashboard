#!/usr/bin/env python3
"""Side product of the trehalose survey: genomes whose Pfam scan, not whose
annotation, is incomplete. BUSCO does not catch these.

Writes data/th_qc_suspect_genomes.tsv for the upstream-follow-ups list in
docs/FUNCTIONAL_ANALYSES.md.
"""
from __future__ import annotations
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"


def main():
    X = pd.read_csv(DATA / "th_genome_matrix.tsv", sep="\t")
    print(f"species reps {len(X)}; median share of proteins with >=1 Pfam hit "
          f"{X.frac_pfam.median():.3f} (q05 {X.frac_pfam.quantile(.05):.3f})")

    print(f"\nreps with no gene_proteins rows at all: {int(X.n_prot.isna().sum())}")
    print(f"reps with 0 of the 1,409 near-universal families: {int((X.n_univ_found == 0).sum())}")

    sus = X[(X.frac_pfam < 0.30) & (X.n_prot >= 500)]
    print(f"\n<30 % of proteins domain-annotated: {len(sus)} reps "
          f"({int((sus.complete_pct >= 90).sum())} of them at BUSCO >= 90)")
    print(sus.CLASS.value_counts().head(8).to_string())

    cols = ["genome", "GENUS", "CLASS", "ord", "SPECIES", "n_prot", "complete_pct", "N50_bp",
            "total_length_bp", "prot_per_mb", "frac_pfam", "n_univ_found", "TPS", "TPP", "HSP70"]
    bad = X[(X.frac_pfam < 0.30) & (X.complete_pct >= 90)].sort_values("frac_pfam")
    bad[cols].to_csv(DATA / "th_qc_suspect_genomes.tsv", sep="\t", index=False)
    print(f"\nwrote data/th_qc_suspect_genomes.tsv ({len(bad)} rows): low domain coverage "
          f"despite BUSCO >= 90")

    # two failure modes, reported separately in the README
    infl = bad[bad.n_prot > 20000]
    scan = bad[bad.n_prot <= 20000]
    print(f"\n  mode (i) partly scanned: {len(scan)} genomes, n_univ_found median "
          f"{scan.n_univ_found.median():.0f}, worst {int(scan.n_univ_found.min())}")
    print(f"  mode (ii) inflated proteome: {len(infl)} genomes "
          f"({', '.join(f'{g} {int(n):,}' for g, n in zip(infl.GENUS, infl.n_prot))}), "
          f"n_univ_found median {infl.n_univ_found.median():.0f}")

    # does the defect touch the trehalose absence calls?
    k = X[X.keep]
    print(f"\nquality-kept TPS-absent genomes: {int((k.TPS == 0).sum())}, of which "
          f"frac_pfam < 0.6: {int(((k.TPS == 0) & (k.frac_pfam < 0.6)).sum())}")
    c = X[(X.TPS == 0) & (X.n_tps.fillna(0) == 0) & (X.qcov >= 0.985)]
    print(f"credible TPS-loss genomes (n={len(c)}): median frac_pfam {c.frac_pfam.median():.3f} "
          f"vs kingdom {X.frac_pfam.median():.3f}; minimum {c.frac_pfam.min():.3f}")
    lic = X[(X.CLASS == "Lecanoromycetes") & X.keep]
    oth = X[X.CLASS.isin(["Dothideomycetes", "Eurotiomycetes", "Leotiomycetes", "Pezizomycetes"]) & X.keep]
    for nm, s in [("Lecanoromycetes", lic), ("Dothio/Eurotio/Leotio/Peziza", oth)]:
        print(f"  {nm:30s} n={len(s):5d} frac_pfam {s.frac_pfam.median():.3f} "
              f"N50 {int(s.N50_bp.median()):>9d} BUSCO {s.complete_pct.median():5.1f} "
              f"TPS-negative {100*(s.TPS == 0).mean():5.1f}%")

    print("\nNote: proteins/Mb is NOT a defect marker on its own. Rusts and ectomycorrhizal "
          "Agaricomycetes genuinely have 100-200 kb-scale gene-poor genomes "
          f"({int(((X.prot_per_mb < 200) & (X.total_length_bp > 5e6)).sum())} reps below 200 "
          "proteins/Mb, most of them real biology). Domain-annotation fraction is the usable "
          "signal because it is measured against the genome's own proteins.")


if __name__ == "__main__":
    main()
