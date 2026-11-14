-- Purpose: rarefied architecture diversity per genome.
-- Inputs: kx2_h1_prot. Output: kx2_h1_arch_metrics.
-- Rarefaction: 2000 proteins per genome chosen by ORDER BY hash(protein_id) (deterministic, no seed needed); genomes with fewer use all.
-- Bigram = unordered pair of adjacent distinct domains. promisc = mean distinct partners over domains with >=1 partner
-- (domains with no partner are excluded). Tandem repeats of the same domain are not counted as a bigram.
-- NOTE: ns counts only sampled proteins WITH a Pfam hit (arch NOT NULL), so per-1000 values in 05 are per 1000 annotated proteins;
-- n_arch/n_bigram are counts within the 2000-protein sample (all 727 genomes with <2000 annotated proteins in the sample = every genome, max ns 1867).
-- Report: distinct architectures/1000, bigrams/1000, promiscuity.
CREATE OR REPLACE TABLE kx2_h1_arch_metrics AS
WITH S AS (SELECT lt, protein_id, arch FROM (SELECT lt, protein_id, arch, row_number() OVER (PARTITION BY lt ORDER BY hash(protein_id)) rn FROM kx2_h1_prot) WHERE rn<=2000),
A AS (SELECT lt, count(*) ns, count(distinct arch) n_arch, count(distinct arch) FILTER (WHERE arch LIKE '%;%') n_arch_multi FROM S WHERE arch IS NOT NULL GROUP BY 1),
U AS (SELECT lt, protein_id, unnest(string_split(arch,';')) AS d, generate_subscripts(string_split(arch,';'),1) AS pos FROM S WHERE arch IS NOT NULL),
B AS (SELECT lt, protein_id, d, lag(d) OVER (PARTITION BY protein_id ORDER BY pos) AS dprev FROM U),
BG AS (SELECT lt, least(d,dprev) a, greatest(d,dprev) b FROM B WHERE dprev IS NOT NULL AND d<>dprev),
BC AS (SELECT lt, count(distinct (a,b)) n_bigram FROM BG GROUP BY 1),
PR AS (SELECT lt, avg(np) promisc, count(*) n_dom_types FROM (SELECT lt, a, count(distinct b) np FROM (SELECT lt,a,b FROM BG UNION ALL SELECT lt,b,a FROM BG) GROUP BY 1,2) GROUP BY 1)
SELECT A.lt AS LOCUSTAG, ns, n_arch, n_arch_multi, coalesce(n_bigram,0) n_bigram, promisc, n_dom_types FROM A LEFT JOIN BC USING(lt) LEFT JOIN PR USING(lt);
