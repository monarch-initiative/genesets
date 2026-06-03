# Expression-Like Query Evals

Disease Ontology sets are useful smoke tests, but they are not realistic query sets for expression-analysis workflows.

For the small curated panel, use `evals/expression20`. This directory is for
larger MyGeneset/MSigDB query snapshots that look more like
differential-expression signatures.

Examples:

```bash
genesets-workflows fetch-mygeneset \
  --query 'GSE*' \
  --source-filter msigdb \
  --limit 2000 \
  --out-dir evals/expression_like/generated/msigdb_gse_2k
```

```bash
genesets-workflows fetch-mygeneset \
  --query 'vs' \
  --source-filter msigdb \
  --limit 5000 \
  --out-dir evals/expression_like/generated/msigdb_vs_5k
```

These query sets can then be used with GO prepared configs by replacing `input.queries` with the generated `queries.gmt`.

The intended scale tests are:

- 2k signatures x GO, single timepoint;
- 10k signatures x GO, single timepoint;
- 2k signatures x 3 GO timepoints x 4 GAF variants;
- larger runs once output filtering and diff summaries are stable.
