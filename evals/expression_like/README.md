# Expression-Like Query Evals

Disease Ontology sets are useful smoke tests, but they are not realistic query sets for expression-analysis workflows.

This directory is for larger MyGeneset/MSigDB query snapshots that look more like differential-expression signatures.

Examples:

```bash
python3 scripts/fetch_mygeneset_query.py \
  --query 'source:msigdb AND gse' \
  --limit 2000 \
  --out-dir evals/expression_like/generated/msigdb_gse_2k
```

```bash
python3 scripts/fetch_mygeneset_query.py \
  --query 'source:msigdb AND name:vs' \
  --limit 5000 \
  --out-dir evals/expression_like/generated/msigdb_vs_5k
```

These query sets can then be used with GO prepared configs by replacing `input.queries` with the generated `queries.gmt`.

The intended scale tests are:

- 2k signatures x GO, single timepoint;
- 10k signatures x GO, single timepoint;
- 2k signatures x 3 GO timepoints x 4 GAF variants;
- larger runs once output filtering and diff summaries are stable.
