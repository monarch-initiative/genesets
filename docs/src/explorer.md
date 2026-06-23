# Web Explorer

The explorer is a local browser for workflow report bundles. It does not run
enrichment itself. It reads `summary.yaml`, query metadata, GMT query genes,
and Parquet result/diff files from one or more existing report directories.

Install or run with the optional explorer dependencies:

```bash
just browser
```

This opens the default set of report bundles, currently:

- the 5,000-set current GOA all-vs-IBA comparison;
- the 5,000-set 2021-vs-2026 GO/GOA comparison;
- the 4,313-set current GOA all-vs-no-`contributes_to` comparison, excluding
  GO-derived query sets.

To open one directly:

```bash
just browser-iba
just browser-go5y
just browser-contributes
```

The equivalent explicit command for the IBA bundle is:

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
```

To open the default browser automatically:

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse --open
```

The server defaults to `http://127.0.0.1:8765`.

## Bundle Inputs

A bundle is any directory with a `summary.yaml` written by `genesets-workflows
go-impact`. The summary points to:

- `queries.gmt`;
- `queries.metadata.json`;
- A-target and B-target result Parquet files;
- threshold-crossing diff Parquet;
- optional term-coverage Parquet files.

If no bundle path is supplied, the explorer looks for
`notebooks/generated/*/summary.yaml`.

## Current Views

The first screen ranks gene sets. The table can be searched and filtered by
source family or stratum, then sorted by criteria such as:

- specific IBA losses;
- all IBA losses;
- specific IBA gains;
- largest p-value delta;
- result rows in either run;
- gene count.

Selecting a gene set opens detail panels:

- `Diffs`: lost, gained, and shared enrichment calls for the configured
  comparison;
- `A targets`: retained enriched targets from the left run;
- `B targets`: retained enriched targets from the right run;
- `Genes`: the query genes from the selected gene set.

The current gene panel shows query genes. To show per-target overlap gene lists,
the workflow layer should add a target-membership artifact such as
`target_gene_membership.parquet` and let the explorer join query genes to target
genes on demand.

## Design Boundary

The web explorer belongs to the workflow layer. It is a reader over report
artifacts and DuckDB queries; it should not add web-specific behavior to the
Rust enrichment kernel.
