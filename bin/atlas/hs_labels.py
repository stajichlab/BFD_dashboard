#!/usr/bin/env python3
"""Build atlas_data/labels/genome_labels_ext.parquet: companion table to genome_labels.parquet
(same LOCUSTAG key, one row per genome) with the columns the size-transition views (HS_*) need.
genome_labels.parquet is not modified, so existing views keep working.

Runs on the login node (small tables only).

Columns added (all genus-level traits unless stated):
  fruitbody_type, hymenium_type, ecm_exploration_type, ecm_lineage, primary_photobiont, decay_type,
  lifestyle_comment                        FungalTraits 1.2 columns not used by genome_labels
  hs_group        organism group for the size-transition views (taxonomy + growth form):
                  yeast, dimorphic, filamentous_other (Pezizomycotina and other non-Agaricomycetes filamentous),
                  Agaricomycetes, lichen, zoosporic, Mucoromycota_sl, Zoopagomycota_sl, Microsporidia, other
  septation       septate_dikarya | coenocytic (Mucoromycota s.l., Zoopagomycota s.l.) | zoosporic | microsporidia | other
  obligate_mode   intracellular | extracellular_obligate | none  (hand list from taxonomy; see OBLIGATE below)
  obligate_lineage  name of the obligate lineage (for per-lineage replicates), else NULL
  micro_size_ok   TRUE when PHYLUM <> Microsporidia, or Microsporidia genome < 10 Mb (PI decision 2026-11-08)
  hs_keep         is_representative AND q_ann AND q_busco_rel AND micro_size_ok (default genome filter of the HS views)
  busco_ge90      busco_complete >= 90 (relaxed family-level filter for Erysiphaceae, Pucciniales, Glomeromycota)
  sister_pair     label of the yeast-origin sister pair the genome belongs to (yeast side or sister side), else NULL
  pair_side       'derived' or 'sister'
"""
import argparse
import os

import duckdb

A = '/bigdata/stajichlab/shared/projects/BFD/BFD_dashboard/atlas_data'

# Obligate biotroph lineages (hand list; taxonomy-based). mode: intracellular vs extracellular_obligate (haustorial or hemocoel).
OBLIGATE = [
    # (lineage name, mode, SQL condition on label columns)
    ('Microsporidia', 'intracellular', "PHYLUM = 'Microsporidia'"),
    ('Rozellomycota', 'intracellular', "PHYLUM = 'Rozellomycota'"),
    ('Aphelidiomycota', 'intracellular', "PHYLUM = 'Aphelidiomycota'"),
    ('Olpidiomycota', 'intracellular', "PHYLUM = 'Olpidiomycota'"),
    ('Pneumocystis', 'intracellular', "GENUS = 'Pneumocystis'"),
    ('Synchytrium', 'intracellular', "GENUS = 'Synchytrium'"),
    ('Pucciniales', 'extracellular_obligate', "ORDER_ = 'Pucciniales'"),
    ('Erysiphaceae', 'extracellular_obligate', "FAMILY = 'Erysiphaceae'"),
    ('Glomeromycota', 'extracellular_obligate', "PHYLUM = 'Glomeromycota'"),
    ('Entomophthorales', 'extracellular_obligate', "ORDER_ = 'Entomophthorales'"),
]

# Yeast-origin sister pairs (derived side, sister side), on label columns.
PAIRS = [
    ('Saccharomycotina vs Pezizomycotina', "yeast_origin = 'Saccharomycotina' AND growth_form = 'yeast'",
     "CLASS IN ('Eurotiomycetes','Dothideomycetes','Sordariomycetes','Leotiomycetes','Pezizomycetes','Orbiliomycetes','Xylonomycetes','Geoglossomycetes') AND growth_form = 'filamentous_mycelium'"),
    ('Taphrinomycotina vs Pezizomycotina', "yeast_origin = 'Taphrinomycotina' AND growth_form = 'yeast'",
     "CLASS IN ('Eurotiomycetes','Dothideomycetes','Sordariomycetes','Leotiomycetes','Pezizomycetes','Orbiliomycetes','Xylonomycetes','Geoglossomycetes') AND growth_form = 'filamentous_mycelium'"),
    ('Tremellomycetes vs Agaricomycetes', "yeast_origin = 'Tremellomycetes' AND growth_form = 'yeast'",
     "CLASS = 'Agaricomycetes' AND growth_form = 'filamentous_mycelium'"),
    ('Pucciniomycotina yeasts vs Pucciniomycetes', "yeast_origin = 'Pucciniomycotina-yeast' AND growth_form = 'yeast'",
     "CLASS = 'Pucciniomycetes'"),
    ('Malasseziomycetes vs Ustilaginomycetes', "CLASS = 'Malasseziomycetes' AND growth_form = 'yeast'",
     "CLASS = 'Ustilaginomycetes'"),
    ('Eurotiomycetes dimorphic vs filamentous', "CLASS = 'Eurotiomycetes' AND growth_form = 'dimorphic_yeast'",
     "CLASS = 'Eurotiomycetes' AND growth_form = 'filamentous_mycelium'"),
]

COENO = "('Mucoromycota','Mortierellomycota','Glomeromycota','Kickxellomycota','Zoopagomycota','Entomophthoromycota','Basidiobolomycota','Calcarisporiellomycota')"
MUCORO = "('Mucoromycota','Mortierellomycota','Glomeromycota')"
ZOOPAGO = "('Kickxellomycota','Zoopagomycota','Entomophthoromycota','Basidiobolomycota')"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--atlas', default=A)
    ap.add_argument('--force', action='store_true')
    a = ap.parse_args()
    out = f'{a.atlas}/labels/genome_labels_ext.parquet'
    if os.path.exists(out) and not a.force:
        print('exists', out); return
    ft = f'{a.atlas}/labels/fungaltraits_1.2.utf8.csv'
    con = duckdb.connect(); con.execute("SET threads=2; SET memory_limit='600MB'")
    con.execute(f"""CREATE TABLE ft AS SELECT GENUS, any_value(Fruitbody_type_template) fruitbody_type, any_value(Hymenium_type_template) hymenium_type,
        any_value(Ectomycorrhiza_exploration_type_template) ecm_exploration_type, any_value(Ectomycorrhiza_lineage_template) ecm_lineage,
        any_value(primary_photobiont) primary_photobiont, any_value(Decay_type_template) decay_type, any_value(Comment_on_lifestyle_template) lifestyle_comment
        FROM read_csv('{ft}', header=true, all_varchar=true) GROUP BY 1""")
    con.execute(f"CREATE TABLE l AS SELECT * FROM '{a.atlas}/labels/genome_labels.parquet'")
    obl_lin = ' '.join(f"WHEN {c} THEN '{n}'" for n, m, c in OBLIGATE)
    obl_mode = ' '.join(f"WHEN {c} THEN '{m}'" for n, m, c in OBLIGATE)
    pair = ' '.join(f"WHEN {d} THEN '{n}' WHEN {s} THEN '{n}'" for n, d, s in PAIRS)
    side = ' '.join(f"WHEN {d} THEN 'derived' WHEN {s} THEN 'sister'" for n, d, s in PAIRS)
    con.execute(f"""CREATE TABLE ext AS SELECT l.LOCUSTAG,
        nullif(ft.fruitbody_type,'') fruitbody_type, nullif(ft.hymenium_type,'') hymenium_type, nullif(ft.ecm_exploration_type,'') ecm_exploration_type,
        nullif(ft.ecm_lineage,'') ecm_lineage, nullif(ft.primary_photobiont,'') primary_photobiont, nullif(ft.decay_type,'') decay_type, nullif(ft.lifestyle_comment,'') lifestyle_comment,
        CASE WHEN PHYLUM = 'Microsporidia' THEN 'Microsporidia'
             WHEN growth_form LIKE 'zoosporic%' OR PHYLUM IN ('Chytridiomycota','Blastocladiomycota','Neocallimastigomycota','Monoblepharomycota','Aphelidiomycota','Rozellomycota','Olpidiomycota','Sanchytriomycota') THEN 'zoosporic'
             WHEN PHYLUM IN {MUCORO} THEN 'Mucoromycota_sl'
             WHEN PHYLUM IN {ZOOPAGO} THEN 'Zoopagomycota_sl'
             WHEN CLASS = 'Agaricomycetes' THEN 'Agaricomycetes'
             WHEN growth_form = 'yeast' THEN 'yeast'
             WHEN growth_form = 'dimorphic_yeast' THEN 'dimorphic'
             WHEN growth_form = 'thallus_photosynthetic' THEN 'lichen'
             WHEN growth_form = 'filamentous_mycelium' THEN 'filamentous_other'
             ELSE 'other' END AS hs_group,
        CASE WHEN PHYLUM = 'Microsporidia' THEN 'microsporidia'
             WHEN PHYLUM IN {COENO} THEN 'coenocytic'
             WHEN PHYLUM IN ('Chytridiomycota','Blastocladiomycota','Neocallimastigomycota','Monoblepharomycota','Aphelidiomycota','Rozellomycota','Olpidiomycota','Sanchytriomycota') THEN 'zoosporic'
             WHEN PHYLUM IN ('Ascomycota','Basidiomycota') THEN 'septate_dikarya' ELSE 'other' END AS septation,
        CASE {obl_mode} ELSE 'none' END AS obligate_mode,
        CASE {obl_lin} ELSE NULL END AS obligate_lineage,
        (PHYLUM IS DISTINCT FROM 'Microsporidia' OR total_length_bp < 10000000) AS micro_size_ok,
        (is_representative AND q_ann AND coalesce(q_busco_rel,false) AND (PHYLUM IS DISTINCT FROM 'Microsporidia' OR total_length_bp < 10000000)) AS hs_keep,
        (busco_complete >= 90) AS busco_ge90,
        CASE {pair} ELSE NULL END AS sister_pair,
        CASE {side} ELSE NULL END AS pair_side
        FROM l LEFT JOIN ft USING (GENUS)""")
    tmp = out + '.tmp'
    con.execute(f"COPY (SELECT * FROM ext ORDER BY LOCUSTAG) TO '{tmp}' (FORMAT PARQUET, COMPRESSION ZSTD)")
    os.replace(tmp, out)
    n = con.execute(f"SELECT count(*) FROM '{out}'").fetchone()[0]
    print('wrote', out, n, 'rows')
    for col in ['hs_group', 'septation', 'obligate_mode', 'obligate_lineage', 'sister_pair']:
        print(col, con.execute(f"SELECT {col}, count(*), count(*) FILTER (WHERE hs_keep) FROM ext e JOIN l USING (LOCUSTAG) GROUP BY 1 ORDER BY 2 DESC").fetchall())


if __name__ == '__main__':
    main()
