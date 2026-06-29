# Schema Reference

The curation schema lives at
`curation/schema/genesets_interpretation.yaml`. It is a LinkML schema for a
single curated interpretation of a non-GO gene set.

## Main Classes

- `GeneSetInterpretation`: one curated gene set interpretation.
- `BiologicalContext`: ontology terms describing what the set is about, such
  as cell type, disease, phenotype, perturbagen, or experimental condition.
- `TermAssociation`: one curated association between the gene set and a GO term.
- `Term`: an ontology term reference with an ID and label.
- `EnrichmentStats`: enrichment values that seeded an association.
- `EvidenceItem`: literature evidence supporting, refuting, or contextualizing
  a term association.

## Important Fields

- `gene_set_id`: stable source identifier, such as `MSIGDB:<SET>`.
- `contexts`: ontology-grounded biological context for the gene set.
- `associations.term`: the GO term being judged.
- `associations.category`: curator role judgment.
- `associations.confidence`: confidence in the biological judgment.
- `associations.specificity`: how specific the term is to the context.
- `associations.recovery_status`: whether current annotations and membership
  make the term recoverable.
- `associations.evidence`: cited evidence with optional validated snippets.

## Validation

Validation has three layers:

- LinkML structural validation checks required fields and enum values.
- `linkml-term-validator` checks ontology IDs and labels for GO, CL, UBERON,
  MONDO, CHEBI, PR, HP, NCBITaxon, EFO, and related prefixes.
- `linkml-reference-validator` checks cited evidence snippets against
  referenced literature when snippets are provided.

Use:

```bash
just curate-validate-schema
just curate-validate
```

## LinkML Browser

A generated LinkML schema browser is useful as reference documentation for the
schema. It should complement, not replace, a curated gene set browser. Most
users will want pages organized by gene set and biological context; schema docs
are organized by classes and slots.
