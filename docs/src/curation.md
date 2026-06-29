# Curated Gene Set Interpretations

The curation corpus stores expert GO interpretations of non-GO gene sets. It is
both a precision/recall fixture for enrichment reports and the source content
for a future curated gene set browser.

Interpretations live in `curation/genesets/*.yaml`. Each file conforms to the
LinkML schema at `curation/schema/genesets_interpretation.yaml` and is validated
with:

- LinkML structural validation;
- `linkml-term-validator` for ontology IDs and labels;
- `linkml-reference-validator` for cited literature snippets.

## What A Curated File Contains

Each curated interpretation records:

- stable source identifier and source collection;
- biological context terms such as cell type, disease, phenotype, perturbagen,
  or experimental condition;
- curated GO term associations;
- curator role category, confidence, and specificity;
- recovery status, distinguishing annotation gaps from gene set membership gaps;
- enrichment stats that seeded an association;
- optional cited evidence and curator notes.

The current corpus includes MSigDB C8 cell-type signatures, KEGG legacy disease
and pathway sets, Reactome and Hallmark sets, and a C7 perturbation contrast.
The manifest at `curation/genesets/manifest.tsv` is the index of curated sets.

## Why This Is Separate From Evals

Eval result tables say which GO terms a run recovered. The curation corpus says
which GO terms are biologically appropriate for a gene set. Keeping those
separate lets reports measure precision and recall without rewriting biological
truth to match the current state of GO annotations.

`category` is the biological judgment. `recovery_status` explains whether the
current annotation and membership state makes that judgment recoverable:

- `annotation_supported`: current annotations recover the term;
- `annotation_gap`: relevant genes are present, but GO annotation is too
  shallow;
- `membership_gap`: the gene set lacks genes needed to support the term.

## Local Commands

Validate the corpus and build the current curation report:

```bash
just curate-validate
just curate-report
```

See [Curation Workflow](curation/workflow.md) for the full workflow and
[Schema Reference](curation/schema.md) for the LinkML model.

## Browser Direction

The curated corpus should become browsable as generated pages or a static site.
The source of truth should remain the validated YAML files, with generated JSON
and HTML/Markdown as read models.

Use LinkML-generated docs for schema reference. Use generated per-gene-set
pages for the corpus itself. See [Static Browser Plan](curation/browser-plan.md).
