# Choosing A Path

Use this page to decide which part of the project to start with.

## I Have Normalized Gene Sets And Want Enrichment

Use the Rust CLI directly:

```bash
genesets-rs matrix \
  --target-sets library.gmt \
  --target-format gmt \
  --queries samples.gmt \
  --query-format gmt \
  --output-format parquet \
  --output results.parquet
```

Start with [Getting Started](getting-started.md), then read the
[Input Model](input-model.md) and [CLI Reference](cli-reference.md).

## I Need To Fetch Sources Or Run A Repeatable Report

Use `genesets-workflows`. The workflow layer owns downloads, source-specific
filtering, GO/GOA or Reactome preparation, metadata, and report summaries:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
```

Start with [Reports And Evals](reports-and-evals.md) and
[Workflow Layer](workflow-layer.md).

## I Want To Inspect Existing Eval Results In A Browser

Use the local web explorer. It reads workflow report bundles and does not run
enrichment itself:

```bash
just browser
```

Start with [Web Explorer](explorer.md) and [Report Bundles](report-bundles.md).

## I Want Curated Gene Set Interpretations

Use the `curation/` corpus. Each curated YAML file records what a non-GO gene
set should mean in GO terms, including curator categories, confidence,
specificity, recovery status, enrichment stats, and evidence.

Start with [Curated Gene Set Interpretations](curation.md), then read the
[Curation Workflow](curation/workflow.md) and [Schema Reference](curation/schema.md).

## I Want A Public Browser For Curated Gene Sets

Build it as static pages or a workflow-layer web view over the curated YAML,
not as a Rust-core feature. The recommended path is described in
[Static Browser Plan](curation/browser-plan.md).
