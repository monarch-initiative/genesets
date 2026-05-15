# Notebook Workflows

The CLI should remain the primary interface, but notebooks are useful for eval review and result exploration.

## Recommended Shape

Use notebooks as consumers of frozen TSV output, not as the primary source of truth.

A typical workflow:

```bash
genesets-rs run evals/disease20/config.yaml \
  --output evals/disease20/generated/results.tsv
```

Then inspect `results.tsv` in a notebook with pandas, Polars, R, or Rust.

## Python Notebook

Python notebooks are the most portable option for collaborators:

- read TSV outputs with pandas or Polars;
- plot p-value distributions;
- compare ranks across ontology releases;
- join results to source metadata.

## Rust Notebook

Rust notebooks are possible with `evcxr_jupyter`:

```bash
cargo install evcxr_jupyter
evcxr_jupyter --install
```

This is useful for demonstrating the library API, but it adds setup friction. The docs should keep CLI-first examples and treat Rust notebooks as optional.

## Repository Policy

Commit small, deterministic notebooks only when they add durable explanation. Large outputs should live under ignored generated directories or external artifacts.
