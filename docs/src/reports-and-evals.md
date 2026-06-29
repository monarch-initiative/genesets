# Reports And Evals

Reports and evals live above the Rust compute engine. They answer questions
that require source-specific preparation, repeated runs, metadata, comparison
logic, and human review.

The main implementation lives in `python/genesets-workflows`. It prepares
inputs, calls `genesets-rs` in coarse batches, writes Parquet and YAML/JSON
artifacts, and then summarizes those artifacts with DuckDB, notebooks, docs, or
the local web explorer.

## What Belongs Here

- fetching public source data from GO, GOA, Reactome, MyGeneSet, or similar;
- applying evidence-code and relation policies;
- preparing normalized term, closure, annotation, query, and background files;
- running configured batches through the Rust CLI;
- comparing ontology releases or annotation variants;
- writing report bundles that can be inspected by SQL, notebooks, or the web UI;
- using curated interpretations as precision/recall fixtures.

## Core Pages

- [Workflow Layer](workflow-layer.md): package boundary and commands.
- [Evals](evals.md): current eval datasets and report findings.
- [Diffing](diffing.md): threshold-crossing comparisons.
- [Post-Processing](post-processing.md): term filtering and interpretation.
- [Tutorials](tutorials.md): runnable shell and notebook workflows.
- [Notebook Workflows](notebooks.md): notebook conventions.
