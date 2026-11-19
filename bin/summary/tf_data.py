"""Loaders and genus aggregation shared by the TF modules (figs_tf.py and later TF statistics).

Tables come from build_tf_tables.py: tf_family_genome (species_prefix, family, n_proteins) and tf_summary_genome
(species_prefix, n_proteins_total, n_tf_proteins, tf_fraction, n_families, shannon), both for every genome.
tf_fraction is NaN for a genome with an empty proteome. shannon is 0 for one family and NaN for none.
"""
from pathlib import Path

import pandas as pd

from suppress import drop_suppressed, load_suppress

FAMILY_MAP = Path(__file__).resolve().parents[2] / "config" / "tf_family_map.tsv"
FAM_PREFIX = "fam_"
PF_PREFIX = "pf_"


def load_tf(data_dir, suppress=None):
    """(tf_family_genome, tf_summary_genome) with suppressed genomes removed. suppress: None (the repo list),
    a path, or a DataFrame with a species_prefix column (an empty one disables)."""
    sup = suppress if isinstance(suppress, pd.DataFrame) else load_suppress(suppress)
    d = Path(data_dir)
    fam = pd.read_parquet(d / "tf_family_genome.parquet")
    summ = pd.read_parquet(d / "tf_summary_genome.parquet")
    return drop_suppressed(fam, sup)[0], drop_suppressed(summ, sup)[0]


def rare_family_labels(accessions, map_path=None):
    """Family labels of the map rows with these Pfam accessions, in the order of `accessions`, without repeats
    (STAT has three accessions and one label). An accession that is not in the map raises KeyError."""
    m = pd.read_csv(map_path or FAMILY_MAP, sep="\t", dtype=str, comment="#", keep_default_na=False)
    lab = dict(zip(m["pfam_acc"], m["family"]))
    out = []
    for a in accessions:
        if a not in lab or not lab[a]:
            raise KeyError(f"Pfam accession {a} is not in {map_path or FAMILY_MAP}")
        if lab[a] not in out:
            out.append(lab[a])
    return out


def genus_table(core_all, fam, summ, busco_min=None, paper_fam=None):
    """One row per genus: the genome with the highest busco_complete (then highest N50_bp, then lowest
    species_prefix) among the genomes with complete TF data (a TF summary row with a proteome and a tf_fraction).
    Genomes with busco_complete < busco_min (or no BUSCO value) are excluded first when busco_min is given.
    Columns: genus, species_prefix, PHYLUM, CLASS, ORDER, busco_complete, n_proteins_total, n_tf_proteins,
    tf_fraction, shannon, and one fam_<family> column per family with the protein count (0 where absent).
    paper_fam (tf_paper_family_genome: species_prefix, paper_family, n_proteins), when given, adds one pf_<paper_family> column
    per paper family with the distinct-protein count (0 where absent).
    Rows with a NULL genus are left out; their number (in core_all after the BUSCO filter) is attrs['n_no_genus']."""
    from summary_ordination import pick_genus_reps   # imported here: it pulls in the ordination libraries
    pool = core_all
    if busco_min is not None:
        pool = pool[pool["busco_complete"] >= busco_min]
    n_no_genus = int(pool["GENUS"].isna().sum())
    ok = summ[(summ["n_proteins_total"] > 0) & summ["tf_fraction"].notna()]["species_prefix"]
    reps = pick_genus_reps(pool[pool["species_prefix"].isin(set(ok))])
    keep = ["species_prefix", "PHYLUM", "CLASS", "ORDER", "busco_complete"]
    t = reps[["GENUS"] + keep].rename(columns={"GENUS": "genus"})
    t = t.merge(summ[["species_prefix", "n_proteins_total", "n_tf_proteins", "tf_fraction", "shannon"]],
                on="species_prefix", how="left")
    wide = fam.pivot_table(index="species_prefix", columns="family", values="n_proteins", aggfunc="sum", fill_value=0)
    wide = wide.reindex(t["species_prefix"]).fillna(0).astype(int)
    wide.columns = [FAM_PREFIX + str(c) for c in wide.columns]
    t = pd.concat([t.reset_index(drop=True), wide.reset_index(drop=True)], axis=1)
    if paper_fam is not None:
        pw = paper_fam.pivot_table(index="species_prefix", columns="paper_family", values="n_proteins", aggfunc="sum",
                                   fill_value=0)
        pw = pw.reindex(t["species_prefix"]).fillna(0).astype(int)
        pw.columns = [PF_PREFIX + str(c) for c in pw.columns]
        t = pd.concat([t.reset_index(drop=True), pw.reset_index(drop=True)], axis=1)
    t = t.sort_values("genus").reset_index(drop=True)
    t.attrs["n_no_genus"] = n_no_genus
    return t
