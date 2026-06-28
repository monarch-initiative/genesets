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
- `genesets/manifest.tsv` - index of curated sets (id, collection, biological context term, series, series_role).
- `genesets/<SET>.yaml` - one interpretation per gene set.
- `cache/` - committed ontology label cache (written by linkml-term-validator).
- `references_cache/` - cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` - sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o genesets/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` - structural + term + reference + obsolescence validation.
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

## Insight (confirmatory vs mechanistic)
`insight` is an optional per-term curator judgment of a GO term's *interpretive
value for this set*, relative to how the set was derived (`collection`,
`context_type`, `direction`). It separates two kinds of correct-and-core term:

- `confirmatory` — entailed by the set's construction or known identity, so it
  recapitulates the design rather than revealing a mechanism (e.g. an
  estradiol-response set enriching for "response to estradiol"; the PRC2-target
  set `BENPORATH_PRC2_TARGETS` enriching for "embryonic morphogenesis").
- `mechanistic` — a non-obvious process that illuminates how/why and is not
  entailed by the set's construction; a genuine enrichment insight.

Free text is allowed (not a closed enum), so curators can record nuance. The
point is to let the eval separate *recall of mechanistic insight* from raw
recall of correct terms — a method that returns only confirmatory terms is
correct but uninformative.

A corpus-wide audit (tightened rule: `mechanistic` only for a specific,
non-obvious process not entailed by the set's construction; all generic
downstream hallmarks — apoptosis, proliferation, generic PI3K/MAPK, disease
OXPHOS-as-known-hallmark, expected immune terms — are `confirmatory`) tags every
insight-bearing term: **760 confirmatory, 86 mechanistic (~90% / 10%)** of 846.
The mechanistic minority concentrates in causal-gene phenotype sets (the
convergent molecular mechanisms behind a phenotype, e.g. RNA splicing / SMN-snRNP
assembly in motor-neuron disease), CRISPR-screen convergences (flavivirus host
factors → ER protein-biogenesis complexes; the *negative*-regulator brake modules
in T-cell and tumor-evasion screens), the down-arm of expression contrasts
(`JISON_SICKLE_CELL_DISEASE_DN`), and a few perturbations. This confirms the
benchmark is predominantly an *annotation/method* test bed; its insight-testing
power is concentrated in a minority of sets.

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

## Pairs and series
Related signatures are linked with two optional fields: `series` (a shared id,
e.g. `SERIES:MICROGLIA_ACTIVATION`) and `series_role` (this set's pole or
position, e.g. `baseline` vs `activated`). One key groups contrasts (DAM vs
homeostatic microglia; A1 neurotoxic vs A2 reparative astrocytes), directional
up/down pairs (`DANG_MYC_TARGETS_UP`/`_DN`), and ordered series — so the eval
can check that opposite poles of one axis resolve to *contrasting* GO
interpretations (a shared term being core-up on one pole and absent or
down on the other). `series_role` is free text because the meaningful poles
differ per series; the `series` id need not be separately defined.

## Evidence
Each `evidence` item carries a `reference` (PMID/DOI), a verbatim `snippet`
(substring-checked against the cited paper by linkml-reference-validator), a
`supports` value (SUPPORT/REFUTE/NEUTRAL), an `explanation`, and an optional
`evidence_source` (HUMAN_CLINICAL / MODEL_ORGANISM / IN_VITRO / COMPUTATIONAL /
OTHER — mirrors dismech for interoperability). A schema rule enforces that any
item with a `snippet` also has a `reference`, so a quote can never bypass
reference validation.

## Term obsolescence
`curate validate` adds a fourth gate: every ontology term id is swept against
its ontology's OAK `obsoletes()` set. This catches obsolete terms the
term-validator misses — the `sqlite:obo:*` builds retain labels for obsolete
classes, so an obsolete id paired with its old label still passes id+label
validation. An obsolete `id:` now fails the build (e.g. `GO:0050663` cytokine
secretion → use `GO:0032635`; `GO:0062023` → `GO:0031012`).
