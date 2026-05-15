# Disease20 vs GO Eval

This eval enriches 20 human Disease Ontology gene sets from MyGeneset.info against official Gene Ontology human annotations.

The prep helper downloads:

- ontology: `https://current.geneontology.org/ontology/go-basic.obo`
- annotations: `https://current.geneontology.org/annotations/goa_human.gaf.gz`

The GO prep script writes ontology-neutral tables for the Rust CLI:

- `terms.tsv`
- `closure.tsv`
- `background_all_goa_symbols.txt`
- one `gene_terms.tsv` and `config.yaml` per annotation-filter variant

`NOT` annotations are always filtered. The starter variants are:

- `all`: all evidence codes, `contributes_to` retained;
- `no_contributes_to`: all evidence codes, `contributes_to` removed;
- `iba`: IBA evidence only;
- `iba_iea`: IBA or IEA evidence only.

The eval uses GAF DB Object Symbol as the gene ID so it can match the MyGeneset.info query GMT symbols. This is convenient for evals but should be documented whenever results are interpreted.

Run the full eval:

```bash
python3 scripts/run_disease20_go_eval.py
```

By default, result TSVs include only rows with Bonferroni-adjusted p-value <= `0.05`. Override with:

```bash
python3 scripts/run_disease20_go_eval.py --max-p-adjust 1.0
```

Run only prep:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/disease20/sets.tsv \
  --out-dir evals/disease20_vs_go/generated

python3 scripts/prepare_go_eval.py \
  --out-dir evals/disease20_vs_go/generated
```

Generated files are ignored by git. The most important outputs are:

- `generated/metadata.yaml`: GO ontology, GAF, closure, background, and filter metadata;
- `generated/run_metadata.yaml`: run commands, p-value cutoff, timing, and per-variant result summaries;
- `generated/<variant>/results.tsv`: enrichment results.
