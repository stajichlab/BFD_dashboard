"""DeepLoc/IDP data audit: coverage, key joins, internal agreement. Read-only on local BFD.duckdb."""
import sys, time
sys.path.insert(0, sys.argv[1] + "/lib")
import bfd_connect
con = bfd_connect.connect(scratch=False, threads=16)
t0 = time.time()
def q(title, sql):
    print("\n##", title, flush=True)
    try:
        r = con.execute(sql)
        cols = [d[0] for d in r.description]
        print("\t".join(cols))
        for row in r.fetchall():
            print("\t".join(str(x) for x in row))
    except Exception as e:
        print("ERROR", e)
    print(f"[{time.time()-t0:.0f}s]", flush=True)

con.execute("CREATE TEMP TABLE gp AS SELECT species_prefix, protein_id, length FROM bfd.gene_proteins")
for t in ["deeploc", "idp_summary", "wolfpsort", "tmhmm", "targetp", "signalp", "predgpi"]:
    q(f"{t}: rows, distinct proteins, distinct species, proteins not in gene_proteins",
      f"""SELECT count(*) AS n_rows, count(DISTINCT protein_id) AS n_prot, count(DISTINCT species_prefix) AS n_species,
          count(*) FILTER (WHERE NOT EXISTS (SELECT 1 FROM gp WHERE gp.protein_id=x.protein_id)) AS n_rows_not_in_gp
          FROM bfd.{t} x""") if t not in ("deeploc", "wolfpsort", "predgpi") else q(
      f"{t}: rows, species, proteins not in gene_proteins (anti-join)",
      f"""SELECT count(*) AS n_rows, count(DISTINCT species_prefix) AS n_species,
          (SELECT count(*) FROM bfd.{t} x ANTI JOIN gp USING (species_prefix, protein_id)) AS n_rows_not_in_gp FROM bfd.{t}""")

q("gene_proteins: rows, species", "SELECT count(*), count(DISTINCT species_prefix) FROM gp")
q("species rows vs gene_proteins species", "SELECT (SELECT count(*) FROM bfd.species) AS species_rows, (SELECT count(DISTINCT species_prefix) FROM gp) AS gp_species")
q("deeploc model values", "SELECT model, count(*) FROM bfd.deeploc GROUP BY 1")
q("deeploc duplicate keys", "SELECT count(*) FROM (SELECT species_prefix, protein_id FROM bfd.deeploc GROUP BY 1,2 HAVING count(*)>1)")
q("deeploc proteins missing from gene_proteins side (gp not in deeploc)",
  "SELECT count(*) FROM gp ANTI JOIN bfd.deeploc USING (species_prefix, protein_id)")
q("idp_summary: proteins in gp lacking a row, by length bin (is absence = no IDP?)",
  """SELECT CASE WHEN g.length<100 THEN '<100' WHEN g.length<300 THEN '100-299' WHEN g.length<1000 THEN '300-999' ELSE '>=1000' END AS lbin,
     count(*) AS n_gp, count(s.protein_id) AS n_with_row, round(100.0*count(s.protein_id)/count(*),1) AS pct_with_row
     FROM gp g LEFT JOIN bfd.idp_summary s USING (species_prefix, protein_id) GROUP BY 1 ORDER BY 1""")
q("idp_summary.length vs gene_proteins.length mismatches",
  "SELECT count(*) AS n, count(*) FILTER (WHERE s.length<>g.length) AS n_mismatch FROM bfd.idp_summary s JOIN gp g USING (species_prefix, protein_id)")
q("idp_summary stats", "SELECT count(*), avg(IDP_fraction), quantile_cont(IDP_fraction,[0.1,0.5,0.9]) FROM bfd.idp_summary")
q("idp segments vs summary: proteins with segments but no summary row",
  "SELECT count(*) FROM (SELECT DISTINCT species_prefix, protein_id FROM bfd.idp) i ANTI JOIN bfd.idp_summary USING (species_prefix, protein_id)")
q("idp segments: sum of lengths vs summary residues mismatches (sample of 5M proteins)",
  """WITH a AS (SELECT species_prefix, protein_id, sum(IDP_length) AS sl FROM bfd.idp WHERE species_prefix < 'F1' GROUP BY 1,2)
     SELECT count(*) AS n, count(*) FILTER (WHERE a.sl<>s.IDP_residues) AS n_mismatch FROM a JOIN bfd.idp_summary s USING (species_prefix, protein_id)""")
q("species with zero rows per tool (vs gene_proteins species)",
  """SELECT 'deeploc' AS tool, count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.deeploc)
     UNION ALL SELECT 'idp_summary', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.idp_summary)
     UNION ALL SELECT 'wolfpsort', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.wolfpsort)
     UNION ALL SELECT 'tmhmm', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.tmhmm)
     UNION ALL SELECT 'predgpi', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.predgpi)
     UNION ALL SELECT 'signalp', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.signalp)
     UNION ALL SELECT 'targetp', count(*) FROM (SELECT DISTINCT species_prefix FROM gp EXCEPT SELECT DISTINCT species_prefix FROM bfd.targetp)""")
q("signalp columns", "DESCRIBE SELECT * FROM bfd.signalp")
q("predgpi feature values", "SELECT feature, source, count(*) FROM bfd.predgpi GROUP BY ALL ORDER BY 3 DESC LIMIT 10")
q("wolfpsort localization values", "SELECT localization, count(*) FROM bfd.wolfpsort GROUP BY 1 ORDER BY 2 DESC")
q("deeploc primary localization (first label)", "SELECT split_part(localizations,'|',1) AS loc, count(*) FROM bfd.deeploc GROUP BY 1 ORDER BY 2 DESC LIMIT 15")

# internal-agreement controls (pre-set before looking): tools agree on well-defined classes
q("CONTROL A: deeploc has_signal_peptide vs wolfpsort extr (confusion)",
  """SELECT d.has_signal_peptide, w.localization='extr' AS wolf_extr, count(*) AS n
     FROM bfd.deeploc d JOIN bfd.wolfpsort w USING (species_prefix, protein_id)
     WHERE d.species_prefix < 'F08' GROUP BY ALL ORDER BY ALL""")
q("CONTROL B: deeploc Transmembrane signal vs tmhmm PredHel>=1 (tmhmm table holds only PredHel>0?)",
  """SELECT d.membrane_types LIKE '%Transmembrane%' AS dl_tm, t.protein_id IS NOT NULL AS has_tmhmm_row, count(*) AS n
     FROM bfd.deeploc d LEFT JOIN bfd.tmhmm t USING (species_prefix, protein_id)
     WHERE d.species_prefix < 'F08' GROUP BY ALL ORDER BY ALL""")
q("CONTROL C: mean IDP fraction by deeploc primary localization (nucleus expected > mitochondrion; sample)",
  """SELECT split_part(d.localizations,'|',1) AS loc, count(*) AS n, round(avg(s.IDP_fraction),3) AS mean_idp
     FROM bfd.deeploc d JOIN bfd.idp_summary s USING (species_prefix, protein_id)
     WHERE d.species_prefix < 'F08' GROUP BY 1 HAVING count(*)>1000 ORDER BY 3 DESC""")
q("CONTROL D: TF proteins (tf_inventory) IDP fraction and nuclear fraction vs all (sample)",
  """SELECT (t.protein_id IS NOT NULL) AS is_tf, count(*) AS n, round(avg(s.IDP_fraction),3) AS mean_idp,
     round(avg((d.localizations LIKE '%Nucleus%')::INT),3) AS frac_nuc
     FROM bfd.deeploc d JOIN bfd.idp_summary s USING (species_prefix, protein_id)
     LEFT JOIN (SELECT DISTINCT species_prefix, protein_id FROM bfd.tf_inventory) t USING (species_prefix, protein_id)
     WHERE d.species_prefix < 'F08' GROUP BY 1""")
