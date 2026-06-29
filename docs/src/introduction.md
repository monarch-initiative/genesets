# Overview

This project is a small monorepo around gene set enrichment, eval exploration,
and curated gene set interpretation.

The main surfaces are:

- `genesets-rs`: a fast Rust compute engine and CLI for ontology-aware and flat
  gene set enrichment.
- `genesets-workflows`: a Python workflow package for source preparation,
  configured evals, reports, metadata, curation helpers, and local tools.
- Web explorer: a local browser for generated eval report bundles.
- Curated interpretation corpus: LinkML-validated YAML interpretations of
  non-GO gene sets, intended both as a gold standard and as browsable content.

The repository should stay together for now. The pieces share fixtures,
schemas, eval definitions, and version history, but each has a clear boundary.

## Core Model

The compute engine is deliberately ontology-neutral. Gene Ontology is the main
near-term target, but the core data model only requires:

- a table of term CURIEs and human-readable names;
- a precomputed child-to-ancestor closure table;
- gene-to-term annotations;
- optional gene names or symbols;
- query sets and background sets.

Flat libraries, such as MSigDB-style collections, are represented as terms with
no closure. User samples are also just flat terms. This keeps the statistical
engine the same for sample-vs-ontology, sample-vs-library, and arbitrary
term-vs-term matrix jobs.

## Boundaries

The Rust core should accept normalized tables and do one job quickly. It should
not know about species, evidence codes, ontology release policy, identifier
mapping, remote services, curation judgments, or web UI state.

Those concerns belong in the workflow and curation layers:

- GO, GOA, Reactome, MyGeneSet, and other source-specific adapters;
- evidence-code filters and release metadata;
- configured report runs and DuckDB summaries;
- local or static web views over generated artifacts;
- LinkML validation and curated interpretation pages.

Start with [Choosing A Path](choosing-a-path.md) if you are not sure which
surface you need.
