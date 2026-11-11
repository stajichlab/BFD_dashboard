# Agent Preferences: Link Conventions

## Preferred Link Targets for Biological Domains and Genes

When generating reports that mention protein domains, gene families, or other biological features, prefer linking to public, curated databases rather than internal descriptions. This ensures:

1. **Consistency** - All users accessing the same links to up-to-date information
2. **Accuracy** - Links point to authoritative, curated resources
3. **Discoverability** - Users can easily explore the referenced entities

### Recommended Link Pattern

For any of the following feature types, use the format:

```
[PF00001](https://www.ebi.ac.uk/interpro/entry/pfam/PF00001/)
```

or equivalently with the full Pfam accession name:

```
[PF00001|Some description](https://www.ebi.ac.uk/interpro/entry/pfam/PF00001/)
```

### Database Preferences

| Feature Type | Preferred Database | Link Format |
|---|---|---|
| **Pfam domains** | Pfam (pfam.xfam.org) | `[PFxxxxx](https://www.ebi.ac.uk/interpro/entry/pfam/PFxxxxx)` |
| **CAZy families** | CAZy (cazy.org) | `[GHxx](https://cazy.org/famGHxx)` or `[GTxx](https://cazy.org/famGTxx)` |
| **MEROPS peptidases** | MEROPS (merops.org) | `[MxxX](https://www.ebi.ac.uk/merops/search/peptidase?name=MxxX)` |
| **UniProt proteins** | UniProt (uniprot.org) | `[GeneID/ProteinName](https://www.uniprot.org/entry/ACC)` |
| **NCBI taxonomy** | NCBI Taxonomy (ncbi.nlm.nih.gov/taxonomy) | `[Taxon name](https://www.ncbi.nlm.nih.gov/taxonomy/?term=)` |

### Example from Current Work

In `docs/analyses/antarctolithica/README.md`, Pfam identifiers were updated from:

```
| PF13812|PPR_3 | 6.5 | 1 | 1 | 4 | 2.91 | 100 |
```

to:

```
| [PF13812](https://www.ebi.ac.uk/interpro/entry/pfam/PF13812) PPR_3 | 6.5 | 1 | 1 | 4 | 2.91 | 100 |
```

This pattern should be applied consistently across all report files in this project whenever domains or gene families are mentioned.

### When NOT to Link

- Brief mentions in running text where a link would disrupt readability
- When the specific database accession is not known or not applicable
- In tables where the link would break the column structure (use inline text instead)

### Implementation Notes

- Always include the `](url)` immediately after the accession identifier
- For table cells, ensure the link doesn't break column parsing (test by rendering the markdown)
- Priority order for links: Pfam > CAZy > MEROPS > UniProt > NCBI
