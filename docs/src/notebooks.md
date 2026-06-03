# Notebook Workflows

The CLI should remain the primary interface, but notebooks are useful for eval review and result exploration.

## Recommended Shape

Use notebooks as consumers of CLI output and report artifacts, not as the
primary source of truth. Scripts and Rust code remain responsible for
reproducible eval generation.

A typical workflow:

```bash
genesets-rs run evals/expression20/config.yaml
```

Then inspect `evals/expression20/generated/results.tsv` in a notebook with
pandas, Polars, DuckDB, or plotting libraries.

For larger temporal reports, prefer a parameterized report config:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
```

The notebook should show the report parameters and generated artifacts, then
query the Parquet outputs.

The committed demonstrator notebooks live in `notebooks/` and are linked from
the [Tutorials](tutorials.md) chapter:

- `01_cli_quickstart.ipynb`;
- `02_expression20_eval.ipynb`;
- `03_go_diff_with_duckdb.ipynb`.

## Python Notebook

Python notebooks are the most portable option for collaborators:

- read TSV outputs with pandas or Polars;
- plot p-value distributions;
- compare ranks across ontology releases;
- join results to source metadata.

The demonstrator notebooks should use `%%bash` for CLI commands and Python for
analysis only. They should not expose Rust APIs.

## Repository Policy

Commit small, deterministic notebooks only when they add durable explanation.
Keep outputs cleared or tiny, and write generated data under
`notebooks/generated/`, ignored eval `generated/` directories, or external
artifacts.
