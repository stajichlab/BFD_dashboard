"""Subphylum is empty (NULL) in both `species` and `samples` for all representative genomes, so the subphylum
level of the summary is derived from the class with this mapping.

Source: the standard fungal classification, written from memory (citation to verify; for example Wijayawardene
et al. / Index Fungorum). The Saccharomycotina and Pezizomycotina class lists also match SACCH and PEZIZO in
bin/atlas/render_views.py. Classes that are not listed fall back to their phylum, so an unlisted class is never
assigned a subphylum silently. Edit this table to correct or extend it; the figures and the group key follow.
"""

CLASS_TO_SUBPHYLUM = {
    # Ascomycota
    **{c: "Saccharomycotina" for c in [
        "Saccharomycetes", "Pichiomycetes", "Dipodascomycetes", "Lipomycetes", "Trigonopsidomycetes",
        "Alloascoideomycetes", "Sporopachydermiomycetes"]},
    **{c: "Pezizomycotina" for c in [
        "Eurotiomycetes", "Dothideomycetes", "Sordariomycetes", "Leotiomycetes", "Pezizomycetes", "Orbiliomycetes",
        "Xylonomycetes", "Geoglossomycetes", "Lecanoromycetes", "Arthoniomycetes", "Lichinomycetes",
        "Laboulbeniomycetes", "Coniocybomycetes", "Collemopsidiomycetes", "Candelariomycetes", "Sareomycetes",
        "Xylobotryomycetes"]},
    **{c: "Taphrinomycotina" for c in [
        "Taphrinomycetes", "Schizosaccharomycetes", "Pneumocystidomycetes", "Pneumocystomycetes", "Neolectomycetes",
        "Archaeorhizomycetes", "Saitoellomycetes"]},
    # Basidiomycota
    **{c: "Agaricomycotina" for c in ["Agaricomycetes", "Tremellomycetes", "Dacrymycetes"]},
    **{c: "Pucciniomycotina" for c in [
        "Pucciniomycetes", "Microbotryomycetes", "Cystobasidiomycetes", "Agaricostilbomycetes", "Atractiellomycetes",
        "Classiculomycetes", "Mixiomycetes", "Tritirachiomycetes"]},
    **{c: "Ustilaginomycotina" for c in ["Ustilaginomycetes", "Exobasidiomycetes", "Malasseziomycetes"]},
}
