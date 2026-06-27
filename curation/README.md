# Curated GO Interpretations of Non-GO Gene Sets

Curated, validated GO-term interpretations of non-GO MSigDB gene sets. The
corpus is deliberately diverse: C8 single-cell cell-type signatures (immune,
epithelial, and hepatic), plus a C2 curated disease pathway and a C7 cytokine
perturbation contrast — exercising the ontology-agnostic `contexts` (cell type,
tissue, disease, chemical perturbagen). Each `genesets/<SET>.yaml` records, per
GO term, a curator judgment (role category, confidence, specificity) and
optional verbatim literature evidence. The data doubles as a precision/recall
gold standard for `genesets-rs`.

## Layout
- `schema/genesets_interpretation.yaml` - LinkML schema (source of truth).
- `conf/` - OAK adapters and reference-validator config.
- `genesets/manifest.tsv` - index of curated sets (id, collection, biological context term).
- `genesets/<SET>.yaml` - one interpretation per gene set.
- `cache/` - committed ontology label cache (written by linkml-term-validator).
- `references_cache/` - cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` - sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o genesets/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` - structural + term + reference validation.
5. `just curate-report` - precision/recall/F1.

## Categories (biology)
`core_process`, `core_component`, `supporting_process`,
`marker_driven_plausible`, `nonspecific`, `false_association`.

The `category` is an assertion about the **biology** and is *not* driven by the
current state of GO annotations. A term that is biologically core stays core
even when no annotation supports it — we never downgrade truth to match
incomplete data.

## Recovery status (gap analysis)
`recovery_status` is **orthogonal** to `category`. It records whether a
biologically-asserted term is actually supported by the gene set's current GO
annotations, so gaps are visible and scoreable without contaminating the
biological judgment:

- `annotation_supported` — genes in the set are annotated to the term; enrichment recovers it.
- `annotation_gap` — relevant genes **are** in the set but GO annotation is too shallow to capture the term. The gap is **GO's** → a GO annotation curation target.
- `membership_gap` — the genes for the term are **not** in the set (the set is core/legacy/incomplete). The gap is the **gene set's**; the term still belongs in a complete set.

This powers two directions: a well-curated set with `core` + `annotation_gap`
finds GO annotation gaps; `core` + `membership_gap` flags gene-set
incompleteness and keeps the eval's recall denominator honest. The report emits
both a biology-complete `recall` and a tool-fair `recall_supportable` (which
excludes membership gaps the set could never recover), plus `core_annotation_gap`
and `core_membership_gap` counts.

See `genesets/KEGG_PARKINSONS_DISEASE.yaml` for the worked example: `neuron
apoptotic process` is an `annotation_gap` (HTRA2/CDK5 present but the death
machinery is annotated only to the generic parent), and `ferroptosis`
(GO:0097707) is a `membership_gap` (real PD biology, but GPX4/ACSL4/SLC7A11 are
not in this legacy KEGG set) — not a `false_association`.

## Sources and the driver-vs-activity axis
Most sets are MSigDB (`MSIGDB:` ids). A second source is **literature-defined**
gene sets (`LIT:` ids, `collection: LIT:DISEASE_ACTIVITY`) taken directly from a
defining paper, where the paper is simultaneously the identity, the membership,
and the evidence. These capture the **biological activity** of a disease or cell
state — the active program a state exhibits (e.g. disease-associated microglia,
neurotoxic A1 reactive astrocytes, the SenMayo senescence/SASP panel) — as
opposed to the mutated **drivers** of a disease (PanelApp/OMIM-style causal
panels, and the `C5:HPO` phenotype gene sets here).

The two need opposite `recovery_status` defaults. In a **driver/causal** panel
the disease's active process is usually an `annotation_gap` (the driver genes
are not annotated to the downstream process they enact). In an **activity**
signature the perturbed processes are `annotation_supported` (the
differentially-expressed genes *are* the readout and carry those GO
annotations), and any driver/risk gene that rides along is at most
`marker_driven_plausible`. A hand-curated activity panel may also legitimately
carry no `nonspecific` housekeeping term at all — e.g. SenMayo contains no
ribosomal-protein genes, so `translation` is not even an enrichment artefact
for it.

## Evidence
Each `evidence` item carries a `reference` (PMID/DOI), a verbatim `snippet`
(substring-checked against the cited paper by linkml-reference-validator), a
`supports` value (SUPPORT/REFUTE/NEUTRAL), an `explanation`, and an optional
`evidence_source` (HUMAN_CLINICAL / MODEL_ORGANISM / IN_VITRO / COMPUTATIONAL /
OTHER — mirrors dismech for interoperability). A schema rule enforces that any
item with a `snippet` also has a `reference`, so a quote can never bypass
reference validation.
