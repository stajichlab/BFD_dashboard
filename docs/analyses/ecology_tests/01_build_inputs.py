#!/usr/bin/env python3
"""Build the unit table and feature matrices for the ecology tests (README.md).

Usage: 01_build_inputs.py OUT_DIR TABLES_DIR LABELS_DB [REPO_ROOT]
  TABLES_DIR  Fungi_BFD_runs/analysis/summary_tables (genome_all, gene_structure_genome, genome_core, cazy_genome, pfam_genome,
              pfam_names, tf_family_genome, tf_summary_genome)
  LABELS_DB   db/BFD_labels.duckdb (v_genome_labels)
Writes units.parquet (one row per genus representative with a FungalTraits lifestyle), features_long.parquet
(species_prefix, feature, rate per 1,000 genes; features pfam:ACC, cazy:FAMILY, tf:FAMILY) and feature_names.parquet.
"""
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = Path(sys.argv[4]) if len(sys.argv) > 4 else HERE.parents[2]
sys.path.insert(0, str(REPO / "bin" / "summary"))

SAPRO = "saprotroph"
MYCORR = {"ectomycorrhizal", "arbuscular_mycorrhizal"}


def coarse_lifestyle(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return None
    if SAPRO in x:
        return "saprotroph"
    if x in MYCORR:
        return "mycorrhizal"
    if x in {"plant_pathogen", "animal_parasite", "mycoparasite", "lichenized"}:
        return x
    return "other"


def main():
    out, tables, labels_db = Path(sys.argv[1]), Path(sys.argv[2]), sys.argv[3]
    out.mkdir(parents=True, exist_ok=True)
    from make_figures import _resolve_suppress
    from suppress import drop_suppressed
    sup = _resolve_suppress(None)
    con = duckdb.connect()
    con.execute(f"ATTACH '{labels_db}' AS lab (READ_ONLY)")
    t = lambda n: f"'{tables}/{n}.parquet'"
    units = con.execute(f"""
        select l.LOCUSTAG as species_prefix, l.PHYLUM, l.CLASS, l."ORDER_" as ORDER_, l.FAMILY, l.GENUS,
               l.primary_lifestyle as lifestyle, l.growth_form,
               a.total_length_bp, a.gc_pct, a.n_genes, a.genes_per_mb, a.coding_fraction, a.median_intron_len_bp,
               g.median_exons_per_gene, g.intron_frac_genes_with_intron, g.median_protein_len_aa,
               c.n_secreted_strict, c.n_merops, c.n_tf_proteins,
               (c.species_prefix is not null) as has_function
        from lab.v_genome_labels l
        join {t('genome_all')} a on a.species_prefix = l.LOCUSTAG
        left join {t('gene_structure_genome')} g on g.species_prefix = l.LOCUSTAG
        left join {t('genome_core')} c on c.species_prefix = l.LOCUSTAG
        where l.is_genus_rep and l.annotation_status = 'annotated' and l.primary_lifestyle is not null""").df()
    if len(sup):
        units, _ = drop_suppressed(units, sup)
    units["lifestyle_coarse"] = units.lifestyle.map(coarse_lifestyle)
    for col in ("total_length_bp", "n_genes", "genes_per_mb", "median_intron_len_bp", "median_protein_len_aa"):
        units[col] = np.log10(units[col].where(units[col] > 0))
    cz = con.execute(f"select species_prefix, cazy_class, sum(n_proteins) n from {t('cazy_genome')} where cazy_class in ('GH','GT','PL','CE','AA','CBM') group by 1,2").df()
    cz = cz.pivot(index="species_prefix", columns="cazy_class", values="n")
    units = units.merge(cz.add_prefix("n_cazy_").reset_index(), on="species_prefix", how="left")
    ng = np.power(10, units.n_genes)
    for cls in ("GH", "GT", "PL", "CE", "AA", "CBM"):
        v = units[f"n_cazy_{cls}"].fillna(0)
        units[f"cazy_{cls}_per1000"] = np.where(units.has_function, 1000 * v / ng, np.nan)
    for src, dst in (("n_secreted_strict", "secretome_per1000"), ("n_merops", "merops_per1000"), ("n_tf_proteins", "tf_per1000")):
        units[dst] = np.where(units.has_function, 1000 * units[src] / ng, np.nan)
    units = units.drop(columns=[c for c in units.columns if c.startswith("n_cazy_")])
    units.to_parquet(out / "units.parquet", index=False)
    ids = set(units.species_prefix)
    ngenes = dict(zip(units.species_prefix, ng))
    feats = []
    pf = pd.read_parquet(tables / "pfam_genome.parquet", columns=["species_prefix", "pfam_acc", "n_proteins"])
    pf = pf[pf.species_prefix.isin(ids)].assign(feature=lambda d: "pfam:" + d.pfam_acc.str.split(".").str[0])
    feats.append(pf.groupby(["species_prefix", "feature"], as_index=False).n_proteins.sum().rename(columns={"n_proteins": "n"}))
    cf = pd.read_parquet(tables / "cazy_genome.parquet", columns=["species_prefix", "cazy_family", "n_proteins"])
    cf = cf[cf.species_prefix.isin(ids) & (cf.cazy_family != "fungi")].assign(feature=lambda d: "cazy:" + d.cazy_family)
    feats.append(cf.groupby(["species_prefix", "feature"], as_index=False).n_proteins.sum().rename(columns={"n_proteins": "n"}))
    tf = pd.read_parquet(tables / "tf_family_genome.parquet", columns=["species_prefix", "family", "n_proteins"])
    tf = tf[tf.species_prefix.isin(ids)].assign(feature=lambda d: "tf:" + d.family.astype(str))
    feats.append(tf.groupby(["species_prefix", "feature"], as_index=False).n_proteins.sum().rename(columns={"n_proteins": "n"}))
    fl = pd.concat(feats, ignore_index=True)
    fl = fl[fl.species_prefix.isin(set(units.species_prefix[units.has_function]))]   # protein tables exist for function-annotated genomes
    fl["rate"] = 1000 * fl.n / fl.species_prefix.map(ngenes)
    fl[["species_prefix", "feature", "rate"]].to_parquet(out / "features_long.parquet", index=False)
    names = pd.read_parquet(tables / "pfam_names.parquet")
    names = pd.DataFrame({"feature": "pfam:" + names.pfam_acc.str.split(".").str[0], "name": names.pfam_id}).drop_duplicates("feature")
    names.to_parquet(out / "feature_names.parquet", index=False)
    print(f"units={len(units)} with_function={int(units.has_function.sum())} feature_rows={len(fl)} features={fl.feature.nunique()}")
    print(units.lifestyle_coarse.value_counts().to_string())


if __name__ == "__main__":
    main()
