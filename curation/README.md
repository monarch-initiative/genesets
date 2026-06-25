# Curated GO Interpretations of Non-GO Gene Sets

Curated, validated GO-term interpretations of non-GO MSigDB gene sets. The
corpus is deliberately diverse: C8 single-cell cell-type signatures (immune,
epithelial, and hepatic), plus a C2 curated disease pathway and a C7 cytokine
perturbation contrast — exercising the ontology-agnostic `contexts` (cell type,
tissue, disease, chemical perturbagen). Each `genesets/<SET>.yaml` records, per
GO term, a curator judgment (role category, confidence, specificity) and
optional verbatim literature evidence. The data doubles as a precision/recall
gold standard for `genesets-rs`.

## Layout
- `schema/genesets_interpretation.yaml` - LinkML schema (source of truth).
- `conf/` - OAK adapters and reference-validator config.
- `genesets/manifest.tsv` - index of curated sets (id, collection, biological context term).
- `genesets/<SET>.yaml` - one interpretation per gene set.
- `cache/` - committed ontology label cache (written by linkml-term-validator).
- `references_cache/` - cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` - sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o genesets/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` - structural + term + reference validation.
5. `just curate-report` - precision/recall/F1.

## Categories
`core_process`, `core_component`, `supporting_process`,
`marker_driven_plausible`, `nonspecific`, `false_association`.
