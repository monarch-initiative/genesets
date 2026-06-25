# Curated GO Interpretations for Non-GO Gene Sets — Design

- **Date:** 2026-06-24
- **Branch:** `feat/curate-genesets`
- **Status:** Design approved, pending spec review

## 1. Problem & Purpose

MSigDB and similar collections contain many gene sets that are *not* derived from
GO — most importantly **C8 single-cell / cell-type signatures**, but also C2
disease/pathway sets, C7 immune sets, and chemical/genetic perturbation
signatures. For such a set (e.g. a dendritic-cell signature) there is a
biologically *correct or reasonable* set of GO terms one would expect a good
enrichment analysis to surface: the processes and components that cell type or
condition actually executes.

Today `docs/src/evals.md` captures this only as informal prose ("`secretory
granule` … fits dendritic/myeloid biology and looks like coverage loss"). We
want to **formalize these judgments as curated, validated data**.

The curated annotations serve two purposes (evaluation-first, but reusable):

1. **Evaluation gold standard** — "what good enrichment should recover." Drives
   precision/recall/F1 for `genesets-rs`, and lets us measure the impact of GO
   version, evidence code (e.g. IBA), and `contributes_to` on *biological
   correctness*, not just raw overlap counts.
2. **Standalone curated resource** — an expert interpretation of each signature,
   with literature evidence, valuable (and potentially citable) independent of
   any one tool.

**Out of scope:** gene sets that are already GO-derived (GO collections). The
schema generalizes across collections via a `collection` field, but the pilot
curates C8.

## 2. Decisions (from brainstorming)

| Decision | Choice |
|---|---|
| Purpose | Both — evaluation-first, structured as a reusable resource |
| Seeding | Enrichment-seeded, then adjudicated (curator confirms/edits, adds missed terms) |
| Assessment categories | 6-value role enum (below) |
| Scoring | Categorical confidence **and** specificity (each: enum + optional float) |
| Pilot scope | Broad C8 batch (dozens of sets) |
| Biological context | Ontology-agnostic (NOT hardwired to CL) — CL/UBERON/MONDO/CHEBI/PR/HP… |
| Home | This repo (`genesets-rs`), reusing dismech's validator tooling pattern |

## 3. Architecture Choice

**Chosen: LinkML-schema'd per-set YAML in `genesets-rs`, validated dismech-style,
with semi-automated drafting.**

Rejected alternatives:
- **Flat TSV/Parquet, no LinkML** — trivially joins eval Parquet, but no nested
  multi-evidence, no term/reference validation, loses the anti-hallucination
  guarantees that motivated citing dismech.
- **Curate inside dismech** — reuses its wired validators, but dismech is
  disease-mechanism scoped; coupling this enrichment eval to an unrelated repo is
  a workaround.

This repo is the right home: it owns the enrichment kernel, the C8 fetch, and the
Parquet/explorer reporting this work plugs into. We borrow only dismech's
*tooling pattern* (LinkML + `linkml-term-validator` + `linkml-reference-validator`
+ OAK label cache).

## 4. Directory Layout

```
curation/
  README.md
  schema/
    genesets_interpretation.yaml          # LinkML schema (source of truth)
  conf/
    oak_config.yaml                        # OAK adapters per prefix (GO, CL, UBERON, MONDO, CHEBI, PR, HP, NCBITaxon)
  .linkml-term-validator.yaml              # term-validator config
  .linkml-reference-validator.yaml         # reference-validator config
  .linkml-reference-validator-sources.yaml # custom reference sources if needed
  cache/                                   # committed ontology label cache (per ontology: terms.csv)
  references_cache/                        # cited-paper cache (gitignored; rebuilt by validator)
  c8/
    manifest.tsv                           # set_id, name, collection, release, n_genes
    HAY_BONE_MARROW_DENDRITIC_CELL.yaml    # one GeneSetInterpretation per file
    ...
```

Conventions (cache layout, validator config file names) mirror dismech so the
tooling drops in unchanged.

## 5. LinkML Schema

`schema/genesets_interpretation.yaml`. Classes:

### `GeneSetInterpretation` (tree root; one per file)
- `gene_set_id` — identifier, e.g. `MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL`
- `gene_set_name` — human-readable
- `collection` — e.g. `C8`, `C2:CP:KEGG_LEGACY`, `C7`
- `msigdb_release` — e.g. `v2025.1.Hs`
- `taxon` — `Term` (NCBITaxon), e.g. `NCBITaxon:9606`
- `direction` — `DirectionEnum` (UP / DN / NA)
- `n_genes` — integer
- `description` — free-text summary of what the signature represents
- `contexts` — multivalued `BiologicalContext` (what the signature is *about*)
- `provenance` — `EnrichmentProvenance` (how the draft was seeded; makes the gold standard reproducible)
- `curation_status` — `CurationStatusEnum` (draft / in_review / reviewed / final)
- `curator`, `curation_date`
- `associations` — multivalued `TermAssociation`

### `BiologicalContext` (ontology-agnostic)
- `term` — `Term {id, label}`, validated by prefix (CL, UBERON, MONDO, DOID, CHEBI, PR, HP, NCBITaxon, EFO…)
- `context_type` — `ContextTypeEnum`: cell_type · anatomical_structure · disease · phenotype · chemical_perturbagen · genetic_perturbation · infection · developmental_stage · cell_line · experimental_condition
- `role_note` — optional string (e.g. "treatment", "vs comparator")

### `TermAssociation`
- `term` — GO `Term {id, label}` (validated)
- `aspect` — `GOAspectEnum`: biological_process (`meaning: GO:0008150`) / cellular_component (`GO:0005575`) / molecular_function (`GO:0003674`)
- `category` — `AssociationCategoryEnum` (see §6)
- `confidence` — `ConfidenceLevelEnum` (high/medium/low); optional `confidence_score` (0–1)
- `specificity` — `SpecificityEnum` (cell_type_specific / lineage_shared / broadly_expressed / ubiquitous); optional `specificity_score` (0–1)
- `seed_source` — `SeedSourceEnum`: `enrichment_recovered` | `curator_added`
  (a `core_*` association with `curator_added` is a **recall miss**)
- `enrichment_stats` — optional `EnrichmentStats` (the seed numbers)
- `evidence` — multivalued `EvidenceItem` (recommended)
- `curator_note` — string

### `EvidenceItem` (mirrors dismech)
- `reference` — PMID/DOI, `implements: linkml:authoritative_reference`
- `reference_title` — string (recommended)
- `supports` — `EvidenceSupportEnum` (SUPPORT / REFUTE / NEUTRAL)
- `snippet` — exact excerpt from the reference, `implements: linkml:excerpt`
  (validated by `linkml-reference-validator`)
- `explanation` — string

### `EnrichmentStats`
- `p_value`, `p_adjust`, `overlap`, `query_size`, `target_size`
- `overlap_genes` — multivalued string

### `EnrichmentProvenance`
- `tool` (e.g. `genesets-rs`), `config_ref` (eval config path/commit)
- `go_version` / ontology snapshot id
- `evidence_codes` — which GO evidence set (all vs IBA-excluded, etc.)
- `background` — background identifier
- `run_date`

### `Term`
- `id` — `uriorcurie`, `identifier: true`
- `label` — string; comment "validated by linkml-term-validator"

### Enums
- `AssociationCategoryEnum` — §6
- `GOAspectEnum`, `DirectionEnum`, `ContextTypeEnum`, `ConfidenceLevelEnum`,
  `SpecificityEnum`, `SeedSourceEnum`, `EvidenceSupportEnum`, `CurationStatusEnum`

## 6. Assessment Category Enum

`AssociationCategoryEnum` (the curator's judgment of each gene-set → GO-term link):

| Value | Meaning |
|---|---|
| `core_process` | Central biological process this context actually executes |
| `core_component` | Characteristic cellular component / structure |
| `supporting_process` | Real, enabling/supporting process (not defining) |
| `marker_driven_plausible` | Plausible given the marker genes but peripheral / generic |
| `nonspecific` | Housekeeping / broadly expressed; not informative for this context |
| `false_association` | Artifactual / biologically wrong for this set |

For eval: **precision** counts `core_*`/`supporting_process` as good and
`false_association`/`nonspecific` as bad among *recovered* terms; **recall**
counts how many `core_*` terms (any `seed_source`) were `enrichment_recovered`.

## 7. Curation Workflow (enrichment-seeded → adjudicated)

New `curate` subcommand group in `python/genesets-workflows` (typer):

1. **`curate draft <set_id>`** — fetch genes → run `genesets-rs enrich` against the
   pinned GO snapshot → emit a draft YAML where every significant term is a
   pre-filled `TermAssociation` (`seed_source: enrichment_recovered`,
   `enrichment_stats` filled, `aspect` auto-derived from GO, `category` left as a
   review placeholder). Fills `provenance`.
2. **LLM-assisted pre-fill (optional, for the dozens-of-sets scale)** — an agent
   proposes `category`/`confidence`/`specificity` plus a candidate PMID + snippet
   per term; output marked `curation_status: draft`. The validators (§8) are the
   guardrail: hallucinated GO IDs and fabricated quotes fail validation.
3. **Human adjudication** — curator confirms/edits categories, adds
   `curator_added` expected-but-missed `core_*` terms (recall), curates evidence,
   sets status to `reviewed`.

## 8. Validation (`just curate-validate`)

Run over the `curation/` tree, added to CI; deps added to `genesets-workflows`:

- `linkml-validate` — structural conformance to the schema.
- `linkml-term-validator validate-data -c conf/oak_config.yaml` — every GO/CL/
  UBERON/MONDO/CHEBI/PR/HP `id`+`label` exists and matches, cached in `cache/`
  (catches AI-hallucinated terms / wrong labels / obsoletes).
- `linkml-reference-validator` — every `snippet` is a real excerpt of its cited
  PMID/DOI, cached in `references_cache/` (catches fabricated literature).

## 9. Eval / Reporting

**`curate report`** joins the gold standard against enrichment output:
- per-set and aggregate **precision** (fraction of recovered terms that are
  core/supporting, not false/nonspecific), **recall** (fraction of curator
  `core_*` terms recovered), **F1**
- stratified by GO version / evidence code (IBA) / `contributes_to`, feeding the
  existing impact analyses and the web explorer
- output: Parquet + summary, consistent with current eval artifacts

## 10. Testing (uv · typer · just · pytest · doctests · CLI tests)

- pytest for the draft generator against a tiny committed fixture gene set
- CLI tests for each `curate` subcommand (typer)
- doctests on schema-helper functions (e.g. precision/recall computation)
- validator smoke over one committed example interpretation
- one hand-curated canonical example checked into `c8/` (e.g. a dendritic-cell
  set) as the reference pattern for curators

## 11. Implementation Phasing (high level; detailed plan via writing-plans)

1. Schema + enums + `conf/oak_config.yaml` + validator config files; one
   hand-authored example YAML; `just curate-validate` green.
2. `curate draft` (mechanical enrichment-seeded draft generation) + tests.
3. Optional LLM-assisted pre-fill step.
4. `curate report` (precision/recall/F1, stratified) + Parquet + explorer hook.
5. Broad C8 batch curation pass.

## 12. Open Questions

- Whether `references_cache/` is committed (dismech keeps a large one) or
  gitignored and rebuilt — default: gitignored, rebuilt by validator.
- Exact float scales for `confidence_score`/`specificity_score` (default: 0–1,
  optional, enum is primary).
