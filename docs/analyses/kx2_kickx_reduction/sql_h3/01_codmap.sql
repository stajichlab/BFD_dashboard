-- 01_codmap.sql
-- Purpose : standard genetic code table (codon -> aa, degeneracy nf, 2-nt box, 3rd base).
-- Inputs  : none.
-- Output  : bfd_scratch.kx2_h3_codmap
-- Notes   : TCAG-ordered standard code string; stops are '*'. nf for '*' is 3 (unused).
-- Report  : used by every codon step (ENC classes, tAI wobble join).
CREATE OR REPLACE TABLE bfd_scratch.kx2_h3_codmap AS
WITH b AS (SELECT unnest(['T','C','A','G']) AS x, unnest([0,1,2,3]) AS i),
c AS (SELECT b1.x||b2.x||b3.x AS codon, b1.i*16+b2.i*4+b3.i AS idx FROM b b1, b b2, b b3),
aa AS (SELECT codon, substr('FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG', idx+1, 1) AS aa FROM c)
SELECT a.codon, a.aa, n.nf, substr(a.codon,1,2) AS box, substr(a.codon,3,1) AS n3
FROM aa a JOIN (SELECT aa, count(*) nf FROM aa GROUP BY 1) n USING (aa);
