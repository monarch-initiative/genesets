# Evals

The eval framework should answer two questions:

1. Are results statistically consistent with trusted implementations on frozen inputs?
2. Is the engine faster on workloads we care about?

## Layers

Unit tests validate local pieces: bitsets, Fisher exact p-values, loaders, and annotation propagation.

CLI integration tests validate the executable surface and TSV output shape.

Eval datasets validate end-to-end workflows. These should be explicit, versioned, and reproducible.

Benchmarks validate speed on synthetic and real workloads.

## Disease20 Starter Eval

The first real-data eval seed is `evals/disease20`. It names 20 human Disease Ontology gene sets from MyGeneset.info. The manifest is intentionally small and hand-reviewed so changes are easy to reason about.

Fetch the current GMT snapshot:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/disease20/sets.tsv \
  --out-dir evals/disease20/generated
```

This writes:

- `queries.gmt`: disease gene sets in GMT format;
- `metadata.json`: source URL, fetch time, gene counts, and skipped IDs.

Generated files are not required for `cargo test`. The goal is to keep normal CI deterministic while making real-data evals one command away.

## Disease20 vs GO

The second eval uses the same 20 MyGeneset.info disease query sets against official GOA human GO annotations.

Run it:

```bash
python3 scripts/run_disease20_go_eval.py
```

This downloads official GO inputs, prepares ontology-neutral tables, runs four annotation-filter variants, and writes companion metadata:

- `all`: all evidence codes, `NOT` filtered, `contributes_to` retained;
- `no_contributes_to`: all evidence codes, `NOT` and `contributes_to` filtered;
- `iba`: IBA only, `NOT` filtered;
- `iba_iea`: IBA or IEA only, `NOT` filtered.

The prep stage writes `metadata.yaml` with source URLs, file digests, GAF header details, closure relation policy, and annotation filter counts. The run stage writes `run_metadata.yaml` with command timings and result summaries.

These filters are intentionally outside the Rust core engine. The core consumes prepared `gene_id, term_id` tables; eval helpers decide how to derive those tables from GAF.

By default, generated eval configs include `max_p_adjust: 0.05`, so result TSVs contain Bonferroni-significant rows. The companion metadata records this cutoff and the per-variant CLI runtime.

## GO Timepoints

The `go_timepoints` eval runs Disease20 against official GO archive releases:

- `2026-03-25`;
- `2021-05-01`;
- `2016-05-01`.

Run it:

```bash
python3 scripts/run_go_timepoints_eval.py
```

The summary metadata contains per-timepoint ontology sizes, GAF metadata, per-variant row counts, run timing, and pairwise diffs over `(query_id, target_id)` result keys.

## Expression-Like Query Sets

Disease Ontology gene sets are convenient but artificial. For scale and realism, use MyGeneset/MSigDB signatures:

```bash
python3 scripts/fetch_mygeneset_query.py \
  --query 'source:msigdb AND gse' \
  --limit 2000 \
  --out-dir evals/expression_like/generated/msigdb_gse_2k
```

MSigDB GSE, `*_UP`, `*_DN`, and `vs` signatures are closer to differential-expression outputs than disease ontology gene sets. These should become the main large-query workloads for timing and GO-version diffing.

## Planned Comparison Evals

The most useful comparison suite will freeze:

- ontology release;
- association release;
- closure relation policy;
- gene universe;
- query gene sets;
- exact correction method;
- expected significant terms from each comparator.

Comparator adapters should run each external tool in containers where possible and emit normalized TSV for comparison.
