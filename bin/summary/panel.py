"""Curated Pfam panel for the functional-domain pages of the Quarto summary.

Verified 2026-11-16 against the data itself: all 57 accessions are present in BFD's pfam table and the
Pfam short name stored there (pfam_id) matches the intended domain (for example PF00176 is now named
SNF2-rel_dom, formerly SNF2_N; PF00154 is RecA_N; PF03443 is AA9, formerly Glyco_hydro_61). The figures
print the name from the data. The *interpretation* of each domain (for example Condensation as an NRPS
signature) is standard literature knowledge, not re-derived here.
"""

PANEL = {
    "Kinases and signalling": {
        "PF00069": "Protein kinase", "PF07714": "Tyr/Ser-Thr kinase", "PF00071": "Ras GTPase",
        "PF00072": "Response regulator", "PF00512": "His kinase A", "PF02518": "HATPase_c",
        "PF00018": "SH3", "PF00169": "PH domain",
    },
    "Transporters": {
        "PF07690": "MFS_1", "PF00005": "ABC transporter ATPase", "PF00083": "Sugar transporter",
    },
    "Cytochrome P450 and redox": {
        "PF00067": "Cytochrome P450", "PF00248": "Aldo/keto reductase", "PF00106": "Short-chain dehydrogenase",
    },
    "Transcription factors": {
        "PF00172": "Zn2Cys6 (Zn_clus)", "PF04082": "Fungal_trans", "PF00096": "zf-C2H2",
        "PF00010": "bHLH", "PF00046": "Homeobox", "PF00170": "bZIP_1", "PF00447": "HSF_DNA-bind",
    },
    "Protein-protein scaffolds": {
        "PF00400": "WD40", "PF00023": "Ankyrin", "PF00515": "TPR_1",
    },
    "Chromatin and DNA metabolism": {
        "PF00176": "SNF2_N", "PF00270": "DEAD helicase", "PF00271": "Helicase_C", "PF00628": "PHD finger",
        "PF00385": "Chromo", "PF00439": "Bromodomain", "PF00154": "RecA",
        "PF02735": "Ku", "PF03730": "Ku_C", "PF03731": "Ku_N",
    },
    "Carbohydrate-active (domain level)": {
        "PF00150": "GH5 cellulase", "PF00704": "GH18 chitinase", "PF00331": "GH10 xylanase",
        "PF00295": "GH28 pectinase", "PF03443": "AA9 LPMO (Glyco_hydro_61)",
    },
    "Peptidases (domain level)": {
        "PF00082": "Subtilisin S8", "PF00026": "Aspartyl A1", "PF00112": "Cysteine C1",
        "PF00246": "Metallocarboxypeptidase M14", "PF00450": "Serine carboxypeptidase S10",
    },
    "Secondary metabolism": {
        "PF00501": "AMP-binding", "PF00109": "Ketoacyl-synt", "PF02801": "Ketoacyl-synt_C",
        "PF00668": "Condensation", "PF00550": "PP-binding", "PF00698": "Acyl_transf_1",
    },
    "Lipid metabolism": {
        "PF01764": "Lipase_3", "PF00561": "Abhydrolase_1",
    },
    "Secreted and effector-like domains": {
        "PF00187": "Chitin_bind_1 (CBM18)", "PF01476": "LysM", "PF00188": "CAP (SCP)", "PF05730": "CFEM",
    },
    "Trehalose and carbon storage": {
        "PF00982": "Glyco_transf_20 (TPS)", "PF02358": "Trehalose_PPase (TPP)", "PF01204": "Trehalase",
        "PF03632": "Glyco_hydro_65m (GH65)",
    },
}

PANEL_ACCESSIONS = [acc for cat in PANEL.values() for acc in cat]
