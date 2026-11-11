import duckdb, os
W=os.path.expanduser('~/bigdata/tmp_kx_scratch/')
con=duckdb.connect()
con.sql("create view g as select * from '/tmp/kx_g_all.parquet'")
con.sql(f"""copy (select *,
 case when PHYLUM='Kickxellomycota' then
        case when GENUS='Coemansia' then 'Kx:Coemansia'
             when ORDER_='Kickxellales' and GENUS in ('Spiromyces','Mycoemilia') then 'Kx:Kickxellales-early'
             when ORDER_='Kickxellales' then 'Kx:Kickxellales-core-other'
             when ORDER_='Harpellales' then 'Kx:Harpellales'
             when ORDER_='Dimargaritales' then 'Kx:Dimargaritales'
             else 'Kx:Ramicandelaber' end
  when PHYLUM='Zoopagomycota' then 'Zoopagales'
  when PHYLUM='Entomophthoromycota' then 'Entomophthorales'
  when PHYLUM='Basidiobolomycota' then 'Basidiobolus'
  when PHYLUM='Mortierellomycota' then 'Mortierellales'
  when ORDER_='Mucorales' then 'Mucorales'
  when ORDER_='Umbelopsidales' then 'Umbelopsidales'
  when ORDER_='Endogonales' then 'Endogonales'
  when PHYLUM='Glomeromycota' then 'Glomeromycota'
  when PHYLUM in ('Chytridiomycota','Monoblepharomycota') then 'Chytrid+Monobleph'
  when PHYLUM='Blastocladiomycota' then 'Blastocladiomycota'
  when PHYLUM='Neocallimastigomycota' then 'Neocallimastigomycota'
  when PHYLUM='Microsporidia' then 'Microsporidia'
  when PHYLUM in ('Rozellomycota','Aphelidiomycota') then 'Rozello+Aphelid'
  when CLASS in ('Saccharomycetes','Pichiomycetes') then 'Sacch:'||CLASS
  when CLASS in ('Dipodascomycetes','Lipomycetes','Trigonopsidomycetes','Alloascoideomycetes','Sporopachydermiomycetes') then 'Sacch:'||CLASS
  when CLASS='Schizosaccharomycetes' then 'Schizosaccharomyces'
  when CLASS='Pneumocystomycetes' then 'Pneumocystis'
  when CLASS='Taphrinomycetes' then 'Taphrinomycetes'
  when CLASS in ('Neolectomycetes','Saitoellomycetes','Novakomycetes') then 'Taphrino-other'
  when CLASS='Malasseziomycetes' then 'Malassezia'
  when CLASS='Ustilaginomycetes' then 'Ustilaginomycetes'
  when CLASS='Exobasidiomycetes' then 'Exobasidiomycetes'
  when CLASS='Wallemiomycetes' then 'Wallemia'
  when CLASS='Tremellomycetes' then 'Tremellomycetes'
  when CLASS='Microbotryomycetes' then 'Microbotryomycetes'
  when CLASS in ('Sordariomycetes','Eurotiomycetes','Dothideomycetes','Leotiomycetes','Pezizomycetes','Orbiliomycetes') then 'Pez:'||CLASS
  when CLASS='Agaricomycetes' then 'Agaricomycetes'
  when CLASS='Pucciniomycetes' then 'Pucciniomycetes'
  else null end as grp
 from g) to '{W}g.parquet'""")
con.sql(f"create view gg as select * from '{W}g.parquet'")
# scan set: species reps for early-diverging + small reduced; genus reps for big Dikarya groups
con.sql(f"""copy (select LOCUSTAG, grp, GENUS, SPECIES, is_representative, is_genus_rep from gg where grp is not null and q_ann and
   ((PHYLUM not in ('Ascomycota','Basidiomycota') and is_representative)
    or (grp in ('Schizosaccharomyces','Pneumocystis','Taphrinomycetes','Taphrino-other','Malassezia','Wallemia','Sacch:Lipomycetes','Sacch:Trigonopsidomycetes','Sacch:Alloascoideomycetes','Sacch:Sporopachydermiomycetes') and is_representative)
    or is_genus_rep)) to '{W}scanset.parquet'""")
print(con.sql(f"select grp,count(*) n,count(distinct GENUS) ng from '{W}scanset.parquet' group by 1 order by 1").df().to_string())
print(con.sql(f"select count(*), min(length(LOCUSTAG)), max(length(LOCUSTAG)) from '{W}scanset.parquet'").fetchall())
