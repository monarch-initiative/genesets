# Workflow Layer

`genesets-rs` is the compute kernel. It should stay small, ontology-neutral,
and fast. The outer convenience layer lives in the Python package under
`python/genesets-workflows`.

The workflow layer owns the parts that are useful but not core enrichment
logic:

- fetching public sources such as GO, GOA, Reactome, and MyGeneSet;
- preparing normalized tables, closures, evidence-code variants, and
  backgrounds;
- running batch analyses through the Rust CLI;
- writing Parquet result artifacts and YAML metadata;
- querying result Parquet with DuckDB;
- computing global report diagnostics such as GO term coverage and terms that
  are scorable but never appear as significant enrichment hits;
- generating report summaries for docs, notebooks, and web views;
- validating curated interpretations and, eventually, generating curated
  corpus browser pages.

## Local Use

Run the packaged workflow CLI from this checkout:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
```

For repeated local use, install it editable:

```bash
python3 -m pip install -e python/genesets-workflows
genesets-workflows doctor
```

The Rust binary still needs to be installed or on `PATH`:

```bash
cargo install --path . --force
```

From outside the repository, use an absolute project path or install the
workflow package into the active environment:

```bash
uv run --project /path/to/genesets-rs/python/genesets-workflows \
  genesets-workflows doctor
```

```bash
python3 -m pip install -e /path/to/genesets-rs/python/genesets-workflows
```

## When To Use Which Layer

Use `genesets-rs` directly when inputs are already normalized and the task is a
single enrichment, matrix enrichment, or result comparison:

```bash
genesets-rs matrix \
  --target-sets reactome.gmt \
  --target-format gmt \
  --queries expression_sets.gmt \
  --query-format gmt \
  --output-format parquet \
  --output results.parquet
```

Use `genesets-workflows` when the task has source-specific prep, multiple Rust
runs, metadata, or report summaries:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
```

The workflow command should usually call Rust once per large batch, not once
per gene set. If a workflow needs thousands of enrichments, it should express
that as a batch plan for the Rust CLI.

Launch the local web explorer over an existing report bundle when the task is
interactive triage rather than batch computation:

```bash
just browser
```

The default browser recipe loads the current GOA all-vs-IBA report, the
2021-vs-2026 GO/GOA temporal report, and the current GOA
all-vs-no-`contributes_to` report. Use `just browser-iba`, `just browser-go5y`,
or `just browser-contributes` to open one analysis directly.

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
```

Use the curation commands when the task is adjudicating or validating curated
GO interpretations:

```bash
uv run --project python/genesets-workflows --extra curation \
  genesets-workflows curate --help
```

## Global Diagnostics

The `go-impact` report also writes term-coverage Parquet files. These join the
prepared GO term table, propagated annotation target sizes, and retained
enrichment results. Each term is classified as:

- `significant`: appears in at least one retained enrichment row;
- `scorable_never_significant`: has propagated annotations in the analysis
  background but never crosses the report threshold for the query collection;
- `unscorable`: has no propagated annotations in the analysis background.

This answers a different question from threshold-crossing diffs: which parts of
the ontology are effectively invisible for this query collection and annotation
variant? The paired coverage file compares those statuses across the two
snapshots or annotation variants.

## Current Commands

GO impact report over two prepared GO snapshots:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
```

Expression signatures against official Reactome as a flat pathway library:

```bash
genesets-workflows reactome-flat
```

Prepare the official Reactome GMT target library:

```bash
genesets-workflows prepare-reactome-flat \
  --out-dir evals/reactome_flat/generated/current
```

Fetch a MyGeneSet query snapshot:

```bash
genesets-workflows fetch-mygeneset \
  --query 'GSE*' \
  --source-filter msigdb \
  --limit 500 \
  --out-dir evals/expression_like/generated/msigdb_gse_500
```

For report-quality eval fixtures, prefer a stratified source plan over the
first `GSE*` hits:

```bash
genesets-workflows fetch-mygeneset-stratified \
  evals/expression_like/msigdb_diverse_5k.yaml
```

The stratified fetcher still emits ordinary `queries.gmt`, `background.txt`,
and `metadata.json` files. The difference is that `metadata.json` records the
quota stratum and search query for each set, which makes later report examples
auditable.

`go-impact` configs can also filter the selected query snapshot with regexes
before running Rust. This is useful for benchmark hygiene, for example excluding
GO-derived MSigDB query sets when the target ontology is GO:

```yaml
query_sets:
  source_dir: expression_like/generated/msigdb_diverse_5k
  limit: 4313
  exclude_regex:
    - "^(GOBP|GOCC|GOMF)_"
```

The report metadata records the selected, considered, and skipped query counts.

Browse one or more generated report bundles in a local web UI:

```bash
genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
```

The old script entry points for these packaged commands remain as
compatibility wrappers around the package modules. Some older eval helpers are
still script-native and should be migrated as the workflow layer expands.

## Interface Boundary

Python should call Rust through batch-oriented commands, not tight subprocess
loops. The desired pattern is:

```text
Python source prep -> normalized tables/YAML plan -> one Rust batch command
-> Parquet outputs -> DuckDB summaries/reports
```

This keeps the subprocess boundary coarse. If a future notebook, web service,
or Python API needs in-process enrichment, a `maturin` binding can wrap the
same Rust core crate without replacing the CLI.

The current choice is CLI-first rather than `maturin`-first because the most
variable code is not the Fisher test. It is fetching, filtering, release
metadata, identifier normalization, report generation, and SQL over Parquet.
Those pieces are faster to evolve in Python. The Rust interface should grow
batch abstractions before Python grows tight bindings.

## Artifact Layout

Workflow outputs should be predictable so notebooks, docs, CI jobs, and web
viewers can consume the same products:

```text
run-dir/
  summary.yaml
  summary.json
  queries.gmt
  queries.metadata.json
  left-results.parquet
  right-results.parquet
  left-vs-right.diff.parquet
  left-vs-right.diff.yaml
```

Use Parquet for durable result tables. Use small YAML/JSON files for metadata,
parameters, source URLs, file hashes, version labels, row counts, and timing.
Use TSV only for small human-inspectable runs.

DuckDB is the default introspection layer over Parquet:

```bash
duckdb -c "SELECT class, count(*) FROM 'left-vs-right.diff.parquet' GROUP BY class"
```

That keeps artifacts immutable and easy to share while preserving SQL
inspection.

## Distribution Plan

The repository should remain a monorepo:

- Rust crate and CLI: crates.io release surface.
- `genesets-workflows`: PyPI release surface.
- Docs, eval configs, and notebooks: same git history and version tags.

The Python package should record the expected `genesets-rs` CLI version and
check it at runtime as the workflow layer matures. Release tags can drive both
crates.io and PyPI publishing.

Short-term distribution:

```bash
cargo install --path .
uv run --project python/genesets-workflows genesets-workflows doctor
```

Longer-term distribution:

```bash
cargo install genesets-rs
python3 -m pip install genesets-workflows
```

If in-process Python calls become necessary, add `maturin` bindings around the
same Rust core crate. The CLI should remain supported because it is the most
transparent interface for batch runs, notebooks, workflow engines, and remote
execution.

## Notebooks And Web Explorer

Notebooks should demonstrate the CLI and analyze generated artifacts. They
should not be the primary workflow engine. A notebook should run a configured
workflow command, then display Parquet-derived tables and plots.

The web explorer follows the same model: select a run directory, read
`summary.yaml`, query Parquet with DuckDB, and render threshold crossings,
timing, top changed targets, enrichment rows, and query genes. No web-specific
behavior enters the Rust enrichment kernel.

The curated gene set browser should follow the same separation. Validate the
YAML corpus, materialize a JSON read model, and generate static pages or a
client-side browser from that model. Do not move curation or browser behavior
into the Rust enrichment core.
