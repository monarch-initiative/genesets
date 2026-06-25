# Curated GO Interpretations

`genesets-rs` curates "correct/reasonable" GO interpretations for non-GO gene
sets (MSigDB C8 cell-type signatures first) as a precision/recall gold standard.

Interpretations live in `curation/genesets/*.yaml`, conform to the LinkML schema at
`curation/schema/genesets_interpretation.yaml`, and are validated by
`linkml-validate`, `linkml-term-validator` (GO/CL/UBERON IDs + labels), and
`linkml-reference-validator` (literature snippets).

Each GO term carries a curator role category (`core_process`, `core_component`,
`supporting_process`, `marker_driven_plausible`, `nonspecific`,
`false_association`), a confidence and a specificity, the seeding enrichment
stats, and optional cited evidence. See `curation/README.md` for the workflow.

Build the report:

    just curate-report
