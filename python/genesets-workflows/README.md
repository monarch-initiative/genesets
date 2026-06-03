# genesets-workflows

`genesets-workflows` is the Python convenience layer around the
`genesets-rs` Rust CLI. It handles repeatable source preparation, batch plans,
report orchestration, metadata, and DuckDB/Parquet summaries. The Rust binary
remains the compute kernel.

Install from this checkout during development:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
```

Or install editable into a Python environment:

```bash
python3 -m pip install -e python/genesets-workflows
genesets-workflows doctor
```

Current commands:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
genesets-workflows reactome-flat
genesets-workflows prepare-reactome-flat --out-dir evals/reactome_flat/generated/current
genesets-workflows fetch-mygeneset --query 'GSE*' --source-filter msigdb --limit 500 --out-dir evals/expression_like/generated/msigdb_gse_500
```

The package intentionally communicates with Rust through batch-oriented CLI
commands that write Parquet. A future `maturin` binding can wrap the same Rust
core if in-process Python calls become necessary.
