# Report Bundles

A report bundle is a directory of immutable artifacts produced by a workflow
run. Bundles are the contract between batch computation, notebook analysis,
docs summaries, and the local web explorer.

The local explorer currently expects `summary.yaml` from `genesets-workflows
go-impact`, plus the files referenced by that summary.

## Minimal Shape

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

The exact filenames may vary by report. `summary.yaml` is the stable entry
point and should record relative paths, source URLs, file digests, parameters,
row counts, timing, and version labels.

## Why Bundles

Bundles keep the interfaces simple:

- Rust writes result tables and comparison outputs.
- Python records source and report metadata.
- DuckDB queries Parquet without loading all rows into memory.
- Notebooks and docs can reproduce summaries from the same artifacts.
- The web explorer can browse a report without rerunning enrichment.

Generated bundles should usually live under an ignored directory such as
`notebooks/generated/` or an explicit external results directory. Small eval
configs and manifests belong in git; large generated report artifacts usually
do not.

## Browser Requirements

For the current explorer, a bundle should provide:

- query genes as GMT;
- query metadata as JSON;
- left and right result Parquet files;
- threshold-crossing diff Parquet;
- optional term-coverage Parquet files.

Future bundle formats should preserve the same design: one small manifest file
points to typed data artifacts, and consumers read the manifest instead of
guessing file names.
