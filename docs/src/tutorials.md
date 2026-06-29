# Tutorials

These tutorials cover compute and report workflows. Each one has a matching
Jupyter notebook under the repository's `notebooks/` directory, but the
commands shown here can also be run directly in a shell.

## Notebook Setup

Install the CLI from the checkout:

```bash
cargo install --path /path/to/genesets-rs --force
```

Install the workflow/notebook environment from the Python package:

```bash
uv run --project python/genesets-workflows --extra notebooks genesets-workflows doctor
```

Generated notebook outputs should be written under `notebooks/generated/`,
which is ignored by git.

## Available Tutorials

- [CLI Quickstart](tutorials/cli-quickstart.md): run the bundled example and
  inspect the enrichment TSV.
- [Expression20 vs GO](tutorials/expression20-eval.md): enrich the small
  expression-derived panel against official GOA human/GO targets.
- [GO Impact Report](tutorials/go-diff-duckdb.md): compare 2021 vs 2026
  GO/GOA over 500 expression-derived gene sets using a parameterized report
  config, then inspect threshold crossings and timings.
