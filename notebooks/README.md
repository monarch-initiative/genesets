# genesets-rs Notebooks

These notebooks demonstrate CLI-first workflows. They use `%%bash` cells for
`genesets-rs` commands and Python only for reading outputs, plotting, and
DuckDB queries.

Recommended setup:

```bash
cargo install --path /path/to/genesets-rs --force
python3 -m pip install jupyter pandas matplotlib duckdb
```

Generated notebook outputs should go under `notebooks/generated/`, which is
ignored by git.

Shared display helpers live in `notebooks/display_utils.py`. Keep table and
metric styling there so tutorial notebooks stay focused on the CLI workflow and
result analysis.

## Notebooks

- `01_cli_quickstart.ipynb`: run the bundled example and inspect the TSV.
- `02_expression20_eval.ipynb`: enrich Expression20 against official GOA
  human/GO targets and summarize the GO terms.
- `03_go_diff_with_duckdb.ipynb`: build a five-year GO impact report over 100
  expression-derived gene sets, including timings and threshold-crossing
  summaries.
