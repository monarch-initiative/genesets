# genesets-rs

`genesets-rs` is a CLI-first, ontology-agnostic gene set enrichment tool. The MVP does one-sided Fisher exact enrichment with optional Bonferroni correction, using dense gene bitsets so repeated term-vs-term and sample-vs-term calculations are mostly popcount work after loading.

The implementation does not hardcode GO. GO is represented as ordinary term IDs, a human-readable term-name table, a precomputed child-to-ancestor closure table, and gene-to-term annotations. If no closure is supplied, the target library is treated as flat.

## Build

```bash
cargo build --release
```

## Documentation

The repository includes an mdBook scaffold for GitHub Pages-style documentation:

```bash
cargo install mdbook
mdbook serve
```

The book source lives in `docs/src`. It covers the CLI, input model, statistics, performance model, eval strategy, ontology prep, and competitive landscape.

## Input Tables

All simple tables may be TSV or CSV. Comment lines beginning with `#` are ignored.

- Terms: `term_curie`, `name`
- Closure: `child_term_curie`, `ancestor_term_curie`; reflexive rows are expected but the loader also adds self-annotation defensively
- Ontology annotations: `gene_curie`, `term_curie`
- Gene names: `gene_curie`, `symbol`
- Background: one `gene_curie` per row

Gene set inputs can be:

- `list`: one gene per row, used as a single set
- `pairwise`: `set_id`, `gene_curie`
- `gene-term`: `gene_curie`, `set_id`
- `gmt`: `set_id`, `description`, `gene...`
- `gmx`: simple GMX, first row set IDs and remaining rows genes
- `gmx-desc`: standard GMX with a description row after the set IDs

## Examples

Single sample against ontology targets:

```bash
cargo run -- enrich \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --sample examples/sample.txt \
  --background examples/background.txt \
  --overlap-genes
```

All ontology terms against themselves:

```bash
cargo run -- matrix \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --queries-from-targets \
  --background examples/background.txt
```

Flat library enrichment, such as MSigDB-style GMT:

```bash
cargo run -- matrix \
  --target-sets path/to/library.gmt \
  --target-format gmt \
  --queries path/to/samples.gmx \
  --query-format gmx \
  --background path/to/background.txt \
  --output enrichments.tsv
```

YAML configuration:

```bash
cargo run -- run examples/enrich.yaml
```

## Output

Output is TSV with:

`query_id`, `query_name`, `target_id`, `target_name`, `overlap`, `query_size`, `target_size`, `background_size`, `p_value`, and adjusted p-value.

With `--overlap-genes`, the CLI also emits `overlap_genes` and `overlap_gene_names`.

## Design Notes

- Terms and query sets are represented as dense `Vec<u64>` bitsets over a shared gene universe.
- Ontology closure is applied during indexing: a gene annotated to a child term is set on every ancestor term.
- Default background is the union of target genes. A supplied background file overrides that and can include genes not annotated to any target.
- Matrix mode parallelizes over queries with Rayon. Output volume can dominate runtime for dense N x N jobs.

## Planned Extensions

- Serialized reusable indexes for repeated GO version comparisons.
- Ranked-list methods alongside Fisher exact.
- Ontologizer-style parent-child, elim, weight, and model-set approaches.
- Convenience wrappers for species/build metadata without baking those concerns into the core model.
- Test-data adapters around mygeneset.info and other public gene set sources.

## Evals

The first real-data eval seed is `evals/disease20`, a manifest of 20 human Disease Ontology gene sets from MyGeneset.info. Fetch a local GMT snapshot with:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/disease20/sets.tsv \
  --out-dir evals/disease20/generated
```

Generated eval files are intentionally separate from `cargo test` so ordinary CI stays deterministic.

The GO eval runner writes significance-focused TSVs by default:

```bash
python3 scripts/run_disease20_go_eval.py --max-p-adjust 0.05
```

Companion YAML metadata records the cutoff, source file digests, prep settings, row counts, and run timing.
