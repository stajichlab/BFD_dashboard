# T3 (xerophile balance): what a xerophily label can be built from

Date: 2026-11-16. Status: **evidence review and design implication; no analysis has been run.** T3 asks whether
TPS / (GSY + GPH) rises with desiccation tolerance, with within-genus contrasts where xerophily varies.

Sources: PubMed (abstracts and metadata retrieved through the NCBI PubMed service, DOIs below) and the project's own report
`sources/report_xerophilic_fungi.md`. Everything below is what the retrieved text says; where a number is quoted it is from an
abstract, not from a table I could read.

## What the literature gives

| Finding | Source |
|---|---|
| Lowest water activity (aw) for cell division in a eukaryote about 0.605, known only for the sugar-tolerant *Xeromyces bisporus* | Stevenson et al. 2015, Environ Microbiol 17:257, PMID 25142751, [doi:10.1111/1462-2920.12598](https://doi.org/10.1111/1462-2920.12598) |
| Established mycelial growth limits about 0.65 aw for *Aspergillus penicilloides* and *X. bisporus*, theoretical limits 0.632 and 0.636; *Wallemia*, *Basipetospora halophila*, *Zygosaccharomyces rouxii*, *Aspergillus* and *Eurotium* named as solute-tolerant eukaryotes (no per-species table in the abstract) | Stevenson et al. 2015, ISME J 9:1333, PMID 25500507, [doi:10.1038/ismej.2014.219](https://doi.org/10.1038/ismej.2014.219) |
| *Aspergillus* sections *Aspergillus* and *Restricti* (subgenus *Aspergillus*) are monophyletic and "comprise xerophilic species that are able to grow on substrates with low water activity"; growth on 0 to 25% NaCl was measured per species | Sklenar et al. 2017, Stud Mycol 88:161, PMID 29158611, [doi:10.1016/j.simyco.2017.09.002](https://doi.org/10.1016/j.simyco.2017.09.002) (full text in PMC5678892, not read here) |
| Section *Restricti* species (*A. penicillioides*, *A. restrictus*, *A. halophilicus*, *A. conicus*, *A. glabripes* and others) are xerophilic isolates from houses, museums and storage rooms | Bastholm et al. 2024, PMID 38857796, [doi:10.1016/j.scitotenv.2024.173880](https://doi.org/10.1016/j.scitotenv.2024.173880); Hagiuda et al. 2022, PMID 35753795, [doi:10.4265/bio.27.65](https://doi.org/10.4265/bio.27.65) |
| Review: xerophilic and xerotolerant means aw < 0.85 (Pitt and Hocking), but "many [*Aspergillus*] can grow up to 0.85", so the review uses 0.80; the ability to grow at lowered aw is a shared character of the xerophilic *Aspergillus* subgenus | Pocsi et al. 2024, Appl Microbiol Biotechnol 108:521, PMID 39560743, [doi:10.1007/s00253-024-13338-5](https://doi.org/10.1007/s00253-024-13338-5) (full text retrieved, but species and section names were stripped from the text, so it cannot be used species by species) |
| Documented low xerotolerance: *Mucor racemosus* "no more xerotolerant than many spoilage bacteria"; *Thamnidium elegans* low xerotolerance; *Penicillium corylophilum* moderately xerotolerant (qualitative, meat isolates) | Lowry and Gill 1984, J Appl Bacteriol 56:193, PMID 6725157, [doi:10.1111/j.1365-2672.1984.tb01339.x](https://doi.org/10.1111/j.1365-2672.1984.tb01339.x) |
| *Z. rouxii* is the main spoilage yeast of a product at aw <= 0.822 (a product aw, not a species minimum) | Iacumin et al. 2021, PMID 34809955, [doi:10.1016/j.fm.2021.103929](https://doi.org/10.1016/j.fm.2021.103929) |

## What it does not give

- **No per-species minimum-aw table** reachable through PubMed abstracts or the full text retrieved. The classic source for
  species-level minima (Pitt and Hocking 1977, J Gen Microbiol 101:35; and the book *Fungi and Food Spoilage*) was not
  retrievable here; the project report cites it, so its per-species statements stay **citation to verify**.
- **No documented non-xerotolerant *Aspergillus* or *Penicillium*.** The 0.85 cutoff in the project report is met by many
  *Aspergillus* species, so a "xerophile" label at 0.85 does not separate species inside the genus, and species with no
  label are unchecked, not non-xerophilic.
- The report's genus-level "Tier 2" flags are, by its own statement, not verified species by species.

## Design implication for T3

1. **Xerophily in *Aspergillus* is clade-structured.** The evidence supports a label at the level of the xerophilic sections
   (sections *Aspergillus* and *Restricti*) versus the rest of the genus, not a per-species scale. That is one or two
   origins, so a within-genus contrast has **one or two independent transitions**: it can describe a difference but cannot
   test it with useful power. Any result must be reported that way.
2. **A species-level, within-genus test is not supportable from these sources.** It would need measured minima for a few
   dozen species of one genus with both tolerant and intolerant members.
3. **What can be built now:** a clade-level label table (xerophilic sections of *Aspergillus*; the Wallemiomycetes, *Xeromyces*,
   *Z. rouxii*, *Debaryomyces hansenii*, *Hortaea werneckii* from the project report) with the PMID of the statement behind
   each row, and a documented set of low-tolerance Mucorales (*Mucor racemosus*, *Thamnidium elegans*) if BFD has them. The
   BFD name-to-section mapping for *Aspergillus* still has to be sourced (the database stores no section).
4. **Quality filter:** budding-yeast genera have no species at BUSCO >= 90 in BFD (for example *S. cerevisiae* best 75.8,
   *Z. rouxii* 76.7), so T3 must use the class-relative quality flag of the atlas, not an absolute cutoff.

## Candidate coverage in BFD (species counts from `v_genome_labels`)

*Aspergillus* 273 species, *Penicillium* 187, *Zygosaccharomyces* 16, *Saccharomyces* 21, *Debaryomyces* 14, *Wallemia* 4,
*Hortaea* 2, *Xeromyces* 1. Not found by name: *A. penicillioides*, *A. restrictus*. Present: *A. glaucus* (1 genome),
*A. wentii* (3), *A. sydowii* (10), *A. niger* (123), *A. fumigatus* (416), *A. nidulans* (19), *A. oryzae* (128), *X. bisporus* (1),
*W. ichthyophaga* (22), *W. mellicola* (27), *Z. rouxii* (40), *D. hansenii* (77), *H. werneckii* (76).

## Files in this folder

Only this file is kept. The earlier name-pattern code (`t3_lib.py`, `t3_design.py`, `test_t3_lib.py`) was removed on 2026-11-21:
its label was a species-name pattern, not an evidence-based xerophily label. It remains in git history at commit c5f80ab.
