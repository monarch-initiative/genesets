# genesets-rs

This repository is a monorepo for gene set enrichment, eval exploration, and
curated gene set interpretation.

- `genesets-rs`: a fast Rust CLI and library for ontology-aware and flat gene
  set enrichment.
- `genesets-workflows`: a Python workflow package for source preparation,
  configured evals, reports, curation helpers, and local browser tools.
- Web explorer: a local, not-yet-deployed UI for browsing generated eval report
  bundles.
- `curation/`: a LinkML-validated corpus of curated GO interpretations for
  non-GO gene sets.

The Rust engine remains ontology-neutral. GO is represented as ordinary term
IDs, a human-readable term-name table, a precomputed child-to-ancestor closure
table, and gene-to-term annotations. If no closure is supplied, the target
library is treated as flat.

## Build The Rust CLI

```bash
cargo build --release
```

## Local Install

For day-to-day use on one machine, install the CLI from this checkout:

```bash
cargo install --path /path/to/genesets-rs --force
```

Cargo installs the binary as `~/.cargo/bin/genesets-rs`. Once that directory is
on your `PATH`, the command works from any directory:

```bash
genesets-rs --help
genesets-rs run /path/to/genesets-rs/examples/enrich.yaml
```

After changing the Rust code, rerun the same `cargo install --path ... --force`
command to refresh the installed binary.

## Documentation

The repository includes an mdBook site for GitHub Pages-style documentation:

```bash
cargo install mdbook
mdbook serve
```

The book source lives in `docs/src`. It is organized around the project
surfaces: compute engine, reports/evals, local web explorer, curated
interpretation corpus, and reference material.

## Workflow Layer

The Rust CLI is the fast enrichment kernel. Repeatable source prep and reports
live in a Python package under `python/genesets-workflows`.

Run the workflow CLI from this checkout:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
```

Or install it editable:

```bash
python3 -m pip install -e python/genesets-workflows
genesets-workflows doctor
```

The package currently exposes:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
genesets-workflows reactome-flat
genesets-workflows prepare-reactome-flat --out-dir evals/reactome_flat/generated/current
genesets-workflows fetch-mygeneset --query 'GSE*' --source-filter msigdb --limit 500 --out-dir evals/expression_like/generated/msigdb_gse_500
genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
genesets-workflows curate --help
```

The older script commands for these workflows remain as compatibility wrappers
around the package; older eval helpers will be migrated incrementally.

The intended boundary is:

```text
Python workflow prep -> normalized tables/YAML plan -> Rust batch command
-> Parquet outputs -> DuckDB/report summaries
```

## Web Explorer

The local explorer browses existing workflow report bundles. It does not run
enrichment itself; it reads `summary.yaml`, query metadata, query genes, and
Parquet result/diff files.

Open the default local bundles:

```bash
just browser
```

Open one bundle explicitly:

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse --open
```

The curated gene set browser should be built as a separate static/generated
view over `curation/genesets/*.yaml`, not as part of the Rust compute engine.

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
genesets-rs enrich \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --sample examples/sample.txt \
  --background examples/background.txt \
  --overlap-genes
```

All ontology terms against themselves:

```bash
genesets-rs matrix \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --queries-from-targets \
  --background examples/background.txt
```

Flat library enrichment, such as MSigDB-style GMT:

```bash
genesets-rs matrix \
  --target-sets path/to/library.gmt \
  --target-format gmt \
  --queries path/to/samples.gmx \
  --query-format gmx \
  --background path/to/background.txt \
  --output enrichments.tsv
```

YAML configuration:

```bash
genesets-rs run examples/enrich.yaml
```

Relative paths inside YAML configs are resolved from the config file's
directory. Relative paths passed directly on the command line are resolved from
the current working directory.

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

The preferred real-data smoke panel is `evals/expression20`, a manifest of 20
human expression-derived MyGeneset.info/MSigDB signatures. Fetch a local GMT
snapshot with:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20/generated
```

Generated eval files are intentionally separate from `cargo test` so ordinary CI stays deterministic.

Run the local overlap smoke:

```bash
genesets-rs run evals/expression20/config.yaml
```

Run the 5-year GO impact report over 500 expression-derived signatures:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
```

That report config runs the core `matrix` and `compare` commands, writes
Parquet outputs plus YAML/JSON metadata, and leaves notebook code to inspect
the generated artifacts rather than rebuild the workflow.

Run the same 500 expression-derived signatures against official Reactome
pathways as a flat GMT target library:

```bash
genesets-workflows reactome-flat
```

The GO eval runner writes significance-focused TSVs by default:

```bash
python3 scripts/run_disease20_go_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20_vs_go/generated \
  --eval-name expression20_vs_go \
  --max-p-adjust 0.05
```

Companion YAML metadata records the cutoff, source file digests, prep settings, row counts, and run timing.

## Curated Gene Set Interpretations

`curation/` holds expert GO-term interpretations of non-GO gene sets, validated
with LinkML term/reference validators and used as a precision/recall gold
standard. Each file records biological context, curated GO associations,
category, confidence, specificity, recovery status, enrichment stats, and
optional literature evidence.

See `curation/README.md` and the mdBook curation section.

    just curate-validate
    just curate-report
