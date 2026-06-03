# Expression20 vs GO Eval

Expression20 vs GO enriches the 20 expression-derived query signatures in
`evals/expression20` against official GOA human annotations and `go-basic.obo`.

Run it:

```bash
python3 scripts/run_disease20_go_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20_vs_go/generated \
  --eval-name expression20_vs_go \
  --description "Twenty expression-derived human MSigDB/MyGeneset.info signatures enriched against official GOA human GO annotations."
```

The helper supports the same annotation variants as the Disease20 GO eval:

- `all`: all evidence codes, `NOT` filtered, `contributes_to` retained;
- `no_contributes_to`: all evidence codes, `NOT` and `contributes_to` filtered;
- `iba`: IBA only, `NOT` filtered;
- `iba_iea`: IBA or IEA only, `NOT` filtered.

Generated files are ignored by git.
