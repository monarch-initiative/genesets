# Expression20 Eval

Expression20 is the preferred small real-query eval panel. It replaces
Disease20 as the main smoke set for ontology-version diff work because the
queries are expression-derived signatures rather than curated disease-gene
knowledge sets.

The manifest contains 20 human MSigDB/MyGeneset.info signatures selected for
workflow breadth:

- in vivo or primary immune cell-type contrasts;
- patient or clinical disease-state contrasts;
- vaccination and infection-response contrasts;
- cytokine-stimulated primary cells;
- drug response and resistance signatures;
- cell-line stimulation/differentiation signatures.

Fetch a local GMT snapshot:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20/generated
```

Run expression-vs-expression overlap enrichment:

```bash
genesets-rs run evals/expression20/config.yaml
```

Run against official GOA human/GO:

```bash
python3 scripts/run_disease20_go_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20_vs_go/generated \
  --eval-name expression20_vs_go \
  --description "Twenty expression-derived human MSigDB/MyGeneset.info signatures enriched against official GOA human GO annotations."
```

The generated files are ignored by git. Commit a generated snapshot only when we
intentionally want a fully frozen benchmark fixture.
