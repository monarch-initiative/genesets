# Compute Engine

The Rust compute engine is the small, fast, ontology-neutral layer. It accepts
normalized tables and gene set files, performs over-representation analysis or
result comparison, and writes durable output formats.

Use this layer when:

- inputs are already normalized;
- you want a single enrichment, matrix run, or result comparison;
- the task should be reproducible from explicit local files;
- workflow-specific concerns such as source download, evidence filtering, and
  report layout are already handled elsewhere.

The core engine does not know about species, GO release policy, evidence codes,
identifier mapping, remote services, curation judgments, or web UI state. Those
belong in the workflow and curation layers.

## Core Pages

- [Getting Started](getting-started.md): install and run the bundled examples.
- [Input Model](input-model.md): normalized table and gene set formats.
- [CLI Reference](cli-reference.md): command surfaces and examples.
- [Statistics](statistics.md): enrichment statistic and correction model.
- [Performance Model](performance.md): bitsets, matrix mode, and output costs.
- [Storage Backends](storage-backends.md): TSV and Parquet output choices.
