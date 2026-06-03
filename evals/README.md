# Evals

This directory contains eval manifests and helper scripts for real-data checks.

Normal `cargo test` does not depend on network access. Evals that fetch public data write to `generated/` directories, which are ignored by git.

## Current Evals

- `expression20`: the preferred small real-query panel, with 20 human MSigDB/MyGeneset.info expression-derived signatures covering cell types, patient disease states, infection/vaccine response, cytokine stimulation, drug response, and cell-line perturbation.
- `expression20_vs_go`: Expression20 enriched against official GOA human GO annotations prepared from `go-basic.obo` and `goa_human.gaf.gz`.
- `disease20`: 20 human Disease Ontology gene sets from MyGeneset.info, retained as an artificial knowledge-curation contrast rather than the core expression workflow eval.
- `disease20_vs_go`: the same 20 query sets enriched against official GOA human GO annotations prepared from `go-basic.obo` and `goa_human.gaf.gz`.
- `go_timepoints`: Disease20 against official GO archive releases from 2026, 2021, and 2016.
- `expression_like`: larger MyGeneset/MSigDB query snapshots, such as GSE or `*_UP`/`*_DN` signatures, for scale tests that look more like differential-expression outputs.
- `reactome_flat`: the 500 expression-like query snapshot enriched against the official Reactome pathway GMT as a flat target library.

## Principles

- Freeze manifests in git.
- Keep generated data out of normal tests.
- Record source URLs and fetch times.
- Prefer simple TSV, GMT, and JSON artifacts.
- Normalize all comparator outputs before comparing result ranks or p-values.
