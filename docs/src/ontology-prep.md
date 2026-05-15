# Ontology Prep

The core engine assumes ontology prep has already happened. That keeps enrichment fast and makes the closure policy explicit.

## Recommended First Path: OAK

OAK, the Ontology Access Kit, is a strong first prep layer because it already handles many biomedical ontology concerns:

- local files and ontology selectors;
- CURIEs, labels, aliases, mappings, and obsolete terms;
- command-line and Python access;
- graph traversal over configured predicates;
- cached downloads for common ontology backends.

OAK's CLI follows the pattern:

```bash
runoak --input HANDLE COMMAND [ARGS]
```

For biological ontologies, OAK supports predicate shortcuts for `is_a` and `part_of`:

```bash
runoak -i sqlite:obo:go ancestors GO:0006915 -p i,p
```

For labels, `fill-table` can populate a table of IDs with labels:

```bash
runoak -i sqlite:obo:go fill-table terms.tsv
```

For this project, the prep target is two TSV files:

```text
terms.tsv:   term_id, name
closure.tsv: child, ancestor
```

The closure should be reflexive and should document which predicates were included.

## When Horned-OWL Makes Sense

Horned-OWL is a Rust OWL library focused on performance and large ontologies. It can parse and manipulate OWL ontologies directly in Rust, and its published docs describe 20x to 40x speedups versus OWL API-based workflows for some large-ontology validation tasks.

It is attractive if we want a fully native prep path, especially for:

- parsing OWL/XML, RDF/XML, Functional Syntax, or Manchester Syntax;
- extracting `SubClassOf` and selected object-property restrictions;
- building relation closure without Python;
- shipping a single Rust binary for prep plus enrichment.

The tradeoff is scope. Correct OWL reasoning is more than graph walking. For MVP prep, OAK or precomputed relation-graph exports are lower risk. A future `genesets-prepare` crate can use Horned-OWL for the fast graph-walk case and clearly label it as such.

## Closure Policy

Every prepared closure should record:

- ontology source and version;
- relation policy, such as `is_a` only or `is_a + part_of`;
- whether inferred axioms were included;
- whether obsolete terms were dropped, migrated, or retained;
- whether the closure is reflexive.

Those details affect enrichment results and must be part of eval metadata.

## GAF Filtering Helpers

GAF filtering belongs in prep, not in the enrichment core. For GO evals, the helper script `scripts/prepare_go_eval.py` implements useful comparison variants:

- always remove `NOT` annotations;
- optionally remove `contributes_to`;
- restrict to IBA;
- restrict to IBA plus IEA.

The script reads GAF qualifier and evidence-code columns directly, writes ordinary `gene_id, term_id` annotation tables, and records the exact filter policy in YAML metadata.
