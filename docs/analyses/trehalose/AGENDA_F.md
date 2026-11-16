# Agenda F: separating gain, loss and shrinkage pressure in reduced genomes

Moved from section F of `docs/FUNCTIONAL_ANALYSES.md` unchanged (plan written 2026-11-10; no compute yet). Items are tracked in `TODO.md`.


**Status:** plan (written 2026-11-10, in the `j11-intron-position` worktree). No new
compute for this section; every item names the table or column it needs.

**Frame.** Analysis A established the phenotype (core Kickxellales: −16.2 Mb,
−4,096 genes, 0.65 against 3.83 introns per transcript, CDS length unchanged,
41 of 2,471 core domains lost, only 4 shared with budding yeasts) and analysis E
(on the unmerged branch) found against multidomain enrichment (28.2 % of proteins
with ≥2 Pfam domains against 33.2 % in Mucorales, δ −0.58) and against yeast-like
loss convergence once the null is stratified by per-family dispensability
(0.93× expected). What neither could do is separate **gene loss** from **annotation
loss**, or **selection on size** from **relaxed constraint**. F1–F12 are the analyses
that can, in the order I would run them.

### Cheap and decisive now (no tree, no new annotation)

**F1. Make every absence call quantitative: a per-genome false-absence model.**
The 1,409 near-universal families in `bfd_scratch.th_univ` give each genome an
empirical absence-propensity score (median 1,407/1,409 found, 5th pct 1,275). Fit
`presence ~ clade + log10(N50) + partial_fraction + busco_complete + qcov + gene_length_bias`
per family (logistic, genus random effect) and report, for each candidate loss, the
predicted false-absence probability next to the observation. Deliverable: a
`pfam_absence_credibility(genome, pfam_acc, p_false_absent)` table. This is agenda
A2's first half and a prerequisite for every other claim in this section.
*Falsifier:* none — it can only lower confidence, never raise it falsely.

**F2. Co-absence modules with a propensity-stratified null (the method that worked in G).**
Obs/exp computed Mantel-Haenszel style across strata of each genome's total missing-family
count, instead of the naive 2×2 Fisher that A used. Apply to the 41 "core losses" of the
Kickxellales, to B's 58 co-absence modules, and to the Ku pair specifically.
*Why it matters:* A's Ku70/Ku80 result is the mechanistic linchpin of the "intron loss was
permissive" story, and A itself notes domains are not independent (Ku + Ku_N lost together).
The stratified null will either firm Ku up or dissolve it.
*Falsifier:* Ku co-absence obs/exp ≤ 3 on the same genome set.

**F3. Are the surviving proteomes *simpler*, or just *smaller*? The dilution test.**
Decompose the proteome-size difference between each reduced lineage and its sister into
(i) families absent, (ii) copies per shared family, (iii) proteins with no Pfam domain.
This is B's un-built G2; E ran a version of it (`kx2_h4_*`: Kickxellales lose 41 % of gene
copies, 59 % of that by contraction; yeasts 44 %, 61 %). F3's addition is to report
**absolute counts of multi-domain proteins vs domain-less proteins**, the framing
HYPOTHESES §P3 wants and E did not: if the multi-domain count is roughly constant while the
domain-less count rises, "simplification" is dilution by orphan and TE-derived models, not
loss of complexity. Needs only `A3/protein_core.parquet`.

**F4. Multidomain, done properly.** E's −5 pp result is confounded both ways: the
multidomain fraction correlates r = +0.72 with BUSCO and −0.70 with partial fraction, and the
comparison group (Mucorales) is post-WGD. Run P3 as designed on
`A3/protein_pfam_arch.parquet` (1.4 GB, chunk like J2b): distinct architectures per 1,000
proteins, distinct unordered domain **pairs** per 1,000, mean distinct domains per Q-model
protein (`n_distinct_domains`, in A3 but dropped from `protein_core`), and a promiscuity
index (mean distinct partner domains per domain, over domains in ≥2 architectures). Controls:
exclude architectures made only of TE domains, stratify by `qcov`, report the
domain-less fraction alongside. *Readout that matters:* if Kickxellales re-pair the same
domain set more often while making fewer proteins, that is genuine architecture-level
compaction; if not, P3 is closed.

**F5. Test shrinkage pressure where shrinkage can be measured within a gene.**
A's `med_cds_ratio ≈ 1.00` is the ortholog-level answer (no CDS shrinkage). The atlas answer
for yeasts is a −2.0 % within-family offset that correlates with genome size *only inside
Saccharomycotina* (ρ +0.44 there, −0.16 in filamentous Pezizomycotina). Run the same V50–V52
machinery (`bin/atlas/h2b_build.py`, `A4/fam_sp_genome`, `A4/genome_family_offset`) on the
Kickxellales and Harpellales, and on Zoopagales as the contrasting route. Two quantities A
and E never separated:
 (a) **within-family CDS length ratio vs genome size within lineage** — is there any slope in
     a lineage whose only remaining lever is protein length?
 (b) **the intergenic half of the story.** E found ~70 % of the −16 Mb is gene number and
     ~22 % tighter intergenic space, and could not measure UTRs because exon length = CDS
     length here. Its convergent-189 bp against divergent-825→561 bp asymmetry is the closest
     thing in the whole Kickxellomycota story to a direct measure of selection on genome size.
     Split every intergenic gap into promoter-side and gene-body-side using
     `gene_intergenic_distances` × `gene_info.strand`, per genus, against Saccharomycotina and
     Mucorales. Promoter spacers compressing more than gene-body spacers = size selection, not
     drift.

### Needs the genus tree (running; PLAN build log 31890419)

**F6. Orthogroups + CAFE/Dollo/Count** (agenda A1) — the only way to get gains, and the only
way to stop calling a contraction a loss. Blocked on `atlas_data/phylo/`.
**F7. PGLS / phylogenetic logistic for Ku presence vs intron density** kingdom-wide (B2),
with independent Ku-loss origins as the effective n.
**F8. Ancestral reconstruction of intron density and of the trehalose/carbon block** — G §4.3
shows the block moves together within class; only a tree says whether it moves together in
time.

### Needs new data

**F9. Repeat classes** (`repeats(genome, class, bp)`, B4). Without it, "genome size reduction"
cannot be split into less DNA per gene versus fewer TEs; `masked_pct` is post-mask and A
already called it weakly informative.
**F10. RNA-seq-supported intron subset** (B3) — A's splice-site and phase table is its least
reliable result, and its spurious-intron worry means the true *Coemansia* count is probably
*below* 0.43.
**F11. tblastn validation queue**, cheapest first: TPS+TPP in Lecanoromycetes and
*Scedosporium* (G §4.5); Ku70/Ku80 in *Coemansia* and Harpellales (A's central claim, at
BUSCO ~71 %); the 41 core Kickxellales losses; then the AA9/tyrosinase/laccase "gains" that E
says are not gains.
**F12. Per-gene codon usage from `gene_CDS`.** E stopped at "no per-gene codon data in `bfd`",
but `gene_CDS` × `gene_transcripts` has everything needed for per-gene GC3, ENC and CAI in one
HPCC job. That turns the tRNA supply-vs-demand result (0.47 against 0.61 in yeasts) into a
test of whether highly expressed gene classes were re-optimized — the actual signature of
translational selection, and the one quantity that distinguishes selection on size from drift.

### What I expect this to leave standing

Reduction in core Kickxellales is mostly **gene number plus intergenic space**, with unchanged
CDS length and no architecture simplification (F3–F5). Yeasts did the same thing with a
different gene set; A and E disagree about convergence only because convergence in *which*
genes is a weaker signal than convergence in *how much*. If that holds, the question stops
being "what did they lose" and becomes "why is the promoter spacer compressible in some
lineages and not others" — which is what F5(b) measures.
