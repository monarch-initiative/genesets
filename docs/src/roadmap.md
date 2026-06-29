# Roadmap

## Short Term

- Add Benjamini-Hochberg FDR.
- Add top-k and p-value-threshold output modes.
- Add serialized prepared indexes.
- Reframe docs around compute, workflows/evals, local explorer, and curated
  interpretation corpus.
- Add eval adapters for GOATOOLS and clusterProfiler.
- Add file digests and ontology metadata joins to `compare` metadata.
- Keep eval runners on Parquet by default.
- Move remaining ad hoc eval/source scripts into `genesets-workflows`.
- Add workflow config schema validation and `genesets-rs` version checks.
- Add a `genesets-workflows curate build-pages` command or equivalent static
  corpus page generator.

## Medium Term

- Add OAK-based ontology prep examples.
- Add native Horned-OWL graph-walk prep as an optional crate or feature.
- Add ranked-list input schema.
- Add ranked-list enrichment statistics.
- Add richer benchmark reports for GitHub Pages.
- Add DuckDB notebooks over Parquet result and diff tables.
- Add hierarchy-aware Reactome prep and compare it with the flat baseline.
- Generate LinkML schema docs and curated gene set instance pages for the docs
  site.
- Add a curated gene set browser over generated corpus JSON.

## Longer Term

- Add topology-aware model families.
- Add richer version-diff reports for ontology release comparisons.
- Add wrapper layers for species, release metadata, and identifier mapping.
- Add visualization notebooks that consume normalized result TSVs.
- Deploy the eval explorer and curated corpus browser from stable generated
  artifacts.
