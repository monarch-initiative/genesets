# Curated GO Interpretations of Non-GO Gene Sets

Curated, validated GO-term interpretations of non-GO MSigDB gene sets
(C8 cell-type signatures first). Each `c8/<SET>.yaml` records, per GO term, a
curator judgment (role category, confidence, specificity), and optional literature
evidence. The data doubles as a
precision/recall gold standard for `genesets-rs`.

## Layout
- `schema/genesets_interpretation.yaml` - LinkML schema (source of truth).
- `conf/` - OAK adapters and reference-validator config.
- `c8/manifest.tsv` - seed list of C8 sets and their biological context terms.
- `c8/<SET>.yaml` - one interpretation per gene set.
- `cache/` - committed ontology label cache (written by linkml-term-validator).
- `references_cache/` - cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` - sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o c8/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` - structural + term + reference validation.
5. `just curate-report` - precision/recall/F1.

## Categories
`core_process`, `core_component`, `supporting_process`,
`marker_driven_plausible`, `nonspecific`, `false_association`.
