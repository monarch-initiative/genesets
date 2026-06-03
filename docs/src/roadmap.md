# Roadmap

## Short Term

- Add Benjamini-Hochberg FDR.
- Add top-k and p-value-threshold output modes.
- Add serialized prepared indexes.
- Add eval adapters for GOATOOLS and clusterProfiler.
- Add file digests and ontology metadata joins to `compare` metadata.
- Keep eval runners on Parquet by default.
- Move remaining ad hoc eval/source scripts into `genesets-workflows`.
- Add workflow config schema validation and `genesets-rs` version checks.

## Medium Term

- Add OAK-based ontology prep examples.
- Add native Horned-OWL graph-walk prep as an optional crate or feature.
- Add ranked-list input schema.
- Add ranked-list enrichment statistics.
- Add richer benchmark reports for GitHub Pages.
- Add DuckDB notebooks over Parquet result and diff tables.
- Add hierarchy-aware Reactome prep and compare it with the flat baseline.

## Longer Term

- Add topology-aware model families.
- Add richer version-diff reports for ontology release comparisons.
- Add wrapper layers for species, release metadata, and identifier mapping.
- Add visualization notebooks that consume normalized result TSVs.
