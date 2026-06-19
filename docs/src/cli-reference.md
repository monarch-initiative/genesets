# CLI Reference

The CLI has four subcommands:

- `enrich`: one query set against many target sets;
- `matrix`: many query sets against many target sets;
- `run`: YAML-configured `enrich` or `matrix`.
- `compare`: threshold-crossing diff between two result tables.

## `enrich`

```bash
genesets-rs enrich \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --sample sample.txt \
  --background background.txt \
  --output results.tsv
```

Use `--target-sets` instead of `--annotations` for flat libraries:

```bash
genesets-rs enrich \
  --target-sets library.gmt \
  --target-format gmt \
  --sample sample.txt
```

Useful options:

- `--sample-format`: `auto`, `list`, `pairwise`, `gene-term`, `gmt`, `gmx`, `gmx-desc`;
- `--sample-set`: select one set from a multi-set sample file;
- `--min-overlap`: suppress rows below an overlap count;
- `--max-p-value`: suppress rows above a raw p-value cutoff;
- `--max-p-adjust`: suppress rows above an adjusted p-value cutoff;
- `--correction`: `bonferroni` or `none`;
- `--output-format`: `tsv`, `parquet`, or `null`; `parquet` requires `--output`,
  and `null` is useful for compute-only profiling;
- `--overlap-genes`: include overlapping gene IDs and names;
- `--threads`: set Rayon worker count.

## `matrix`

```bash
genesets-rs matrix \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --queries queries.gmx \
  --query-format gmx \
  --background background.txt
```

Use targets as queries for term-vs-term runs:

```bash
genesets-rs matrix \
  --annotations gene_terms.tsv \
  --closure closure.tsv \
  --queries-from-targets
```

## `run`

```yaml
mode: enrich
ontology:
  terms: terms.tsv
  closure: closure.tsv
  annotations: gene_terms.tsv
input:
  sample: sample.txt
  sample_format: list
  sample_name: sample
background:
  file: background.txt
overlap_genes: true
max_p_adjust: 0.05
output_format: tsv
```

Run it:

```bash
genesets-rs run examples/enrich.yaml
```

Relative paths in YAML are resolved from the config file's directory. Relative
paths passed directly as CLI arguments are resolved from the current working
directory.

For mass evals, write Parquet and inspect it with DuckDB:

```bash
genesets-rs matrix ... --output-format parquet --output results.parquet
duckdb -c "SELECT * FROM 'results.parquet' WHERE p_adjust_bonferroni <= 0.05"
```

## `compare`

Compare two enrichment result tables by `(query_id, target_id)` and classify
adjusted p-value threshold crossings:

```bash
genesets-rs compare \
  --left go-2021.parquet \
  --right go-2026.parquet \
  --p-adjust-cutoff 0.05 \
  --output-format parquet \
  --output go-2021-vs-2026.diff.parquet \
  --metadata-output go-2021-vs-2026.diff.yaml
```

Input formats are inferred from `.tsv`, `.txt`, `.parquet`, or `.pq`, or can be
set explicitly with `--left-format` and `--right-format`. TSV output goes to
stdout by default. Parquet output requires `--output`.

The default output includes:

- `lost_significant`: significant on the left, not significant on the right;
- `gained_significant`: not significant on the left, significant on the right;
- `shared_significant`: significant on both sides.

Use `--crossings-only` to emit only gained/lost rows.

## Workflow CLI

The separate `genesets-workflows` command is the convenience layer for
configured source prep and reports. It calls the Rust CLI for batch compute and
then writes Parquet plus metadata:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
genesets-workflows reactome-flat
```

For interactive review of an existing report bundle, use the optional web
explorer dependencies:

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
```

Use `genesets-rs` for normalized single jobs. Use `genesets-workflows` when the
task needs downloads, evidence filters, release metadata, multiple Rust runs,
DuckDB summaries, notebook/report artifacts, or browser-based result triage.
