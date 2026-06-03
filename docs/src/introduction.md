# Introduction

`genesets-rs` is a fast, CLI-first enrichment engine for ontology-aware and flat gene set analysis.

The project is deliberately ontology-neutral. Gene Ontology is the main near-term use case, but the core data model only requires:

- a table of term CURIEs and human-readable names,
- a precomputed child-to-ancestor closure table,
- gene-to-term annotations,
- optional gene names or symbols,
- query sets and background sets.

Flat libraries, such as MSigDB-style collections, are represented as terms with no closure. User samples are also just flat terms. This keeps the statistical engine the same for sample-vs-ontology, sample-vs-library, and arbitrary term-vs-term matrix jobs.

## What This Is Optimized For

The initial design is aimed at repeated over-representation analysis:

- compare thousands of disease, pathway, or experiment-derived gene sets against one ontology;
- compare two ontology releases by running the same query set against both;
- run ontology terms against themselves to inspect term overlap structure;
- cache or reuse prepared annotation tables from external ontology tooling.

The MVP statistic is one-sided Fisher exact enrichment with optional Bonferroni correction. The implementation leaves room for ranked-list statistics and ontology-aware model-set methods, but the first requirement is a very fast and auditable baseline.

## Project Layers

The project has two user-facing layers:

- `genesets-rs`: the Rust CLI and library code for loading normalized tables,
  building bitset indexes, computing enrichment, comparing result tables, and
  writing TSV or Parquet.
- `genesets-workflows`: the Python workflow package for repeatable source
  preparation, configured reports, metadata, DuckDB summaries, notebooks, and
  future web-facing reports.

The boundary is intentional. GO, GOA, Reactome, MyGeneSet, evidence-code
filters, release metadata, and report layouts change faster than the enrichment
kernel. They belong in the workflow layer. The Rust core should continue to
accept normalized ontology-neutral tables and do one job quickly.

## Non-Goals For The Core Engine

The core engine does not know about species, evidence codes, ontology release policy, identifier mapping, or remote services. Those concerns belong in prep and wrapper layers. This separation keeps the fast path small and makes it easier to build reproducible evals.
