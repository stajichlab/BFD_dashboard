#!/usr/bin/env python3
"""Step 4b: turn one miniprot GFF into per-locus protein models for a genome and family.

Usage: 04b_translate_models.py <genome.fa(.gz)> <locustag> <family> <miniprot.gff> <out_prefix>
Keeps the best-scoring model per locus (fam_lib.greedy_loci), translates it from the genome with the standard genetic code
(fam_lib.model_protein, tested in tests/), and writes <out_prefix>.faa (models passing the filters) and <out_prefix>.tsv (all loci).
Filters, fixed before any result was seen:
  - no frameshift (a change of reading frame between CDS blocks, found from the GFF3 phases) and no internal stop codon (such a model is a
    pseudogene or a mis-spliced or frameshifted alignment; it stays in the table with its status but is not used)
  - length >= MIN_LEN[family] aa (a fraction of the family core domain: ligase 300, RecA/Rad51 200, Pol X 200)
Header of kept models: >MP|<locustag>|<scaffold>:<start>-<end>:<strand>|<best_reference_clade>   (the reference clade is informational
only; classification uses the tree).
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import fam_lib as F

MIN_LEN = {'lig': 300, 'reca': 200, 'polx': 200}
MAX_INTRON = 5000          # models with a longer intron are flagged 'long_intron' and not used (miniprot is also run with -G 5k)
genome_path, loc, fam, gff, prefix = sys.argv[1:6]
genome = F.read_fasta(genome_path)
models = F.parse_miniprot_gff(open(gff))
loci = F.greedy_loci(models)
n_kept = 0
with open(prefix + '.faa', 'w') as fa, open(prefix + '.tsv', 'w') as tsv:
    tsv.write('family\tlocustag\tmodel\tscaffold\tstrand\tstart\tend\tscore\tidentity\tbest_ref\tprot_len\tinternal_stops\tstatus\tmax_intron\tframe_breaks\n')
    for m in loci:
        prot, stops, ok, breaks = F.model_protein(m, genome)
        mi = F.max_intron(m)
        if not ok:
            status = 'no_sequence'
        elif mi > MAX_INTRON:
            status = 'long_intron'
        elif breaks > 0:
            status = 'frameshift'              # the alignment changed reading frame between blocks (miniprot 'Frameshift=true')
        elif stops > 0:
            status = 'internal_stop'
        elif len(prot) < MIN_LEN[fam]:
            status = 'too_short'
        else:
            status = 'kept'
        tsv.write(f"{fam}\t{loc}\t{m['id']}\t{m['scaf']}\t{m['strand']}\t{m['start']}\t{m['end']}\t{m['score']}\t{m['identity']}\t{m['target']}\t{len(prot)}\t{stops}\t{status}\t{mi}\t{breaks}\n")
        if status == 'kept':
            clade = m['target'].split('|')[0]              # reference clade, or 'ANN' when the best query was an annotated BFD protein
            fa.write(f">MP|{loc}|{m['scaf']}:{m['start']}-{m['end']}:{m['strand']}|{clade}\n{prot}\n")
            n_kept += 1
print(f'{loc} {fam}: {len(models)} raw models, {len(loci)} loci, {n_kept} kept')
