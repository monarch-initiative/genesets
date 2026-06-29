# Project Map

This repository is a monorepo for four related surfaces:

| Surface | Main paths | Audience | Responsibility |
| --- | --- | --- | --- |
| Rust compute engine | `src/`, `Cargo.toml` | CLI users, workflow runners, library contributors | Load normalized gene set tables, build bitset indexes, compute enrichment and comparison results, write TSV or Parquet. |
| Python workflows | `python/genesets-workflows/` | Report authors, evaluators, data-prep users | Fetch public sources, prepare normalized inputs, run configured batches, summarize Parquet outputs, and launch local tools. |
| Eval explorer | `python/genesets-workflows/src/genesets_workflows/explorer/` | People triaging report outputs | Browse existing workflow report bundles in a local web UI. |
| Curated interpretation corpus | `curation/` | Curators, evaluators, browser users | Store LinkML-validated GO interpretations of non-GO gene sets, with curator judgments, recovery status, and literature evidence. |

The Rust crate should stay focused on the compute model. GO, GOA, Reactome,
MyGeneSet, evidence-code policies, report layouts, web views, and curation
workflow are deliberately outside the Rust core because they change faster and
are easier to audit as files and workflow artifacts.

## Repository Layout

```text
src/
  Rust library and CLI implementation.
python/genesets-workflows/
  Python package for source prep, reports, curation helpers, and the local explorer.
evals/
  Small and parameterized eval definitions.
curation/
  LinkML schema, validator config, manifest, and curated interpretation YAML.
docs/src/
  mdBook documentation.
notebooks/
  Analysis notebooks that consume generated workflow artifacts.
```

## Release Boundaries

Keep one repository for now, with separate release surfaces:

- `genesets-rs`: crates.io or `cargo install` release surface.
- `genesets-workflows`: PyPI or editable Python package release surface.
- curated interpretation corpus: versioned data release tied to repository tags.
- docs and static generated browser pages: GitHub Pages or another static host.

Split the repository only if one of those surfaces needs an independent
governance, release cadence, or storage policy that makes the monorepo painful.
That is not true yet.
