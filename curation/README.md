# Curated Gene Set Interpretations

Curated, validated GO-term interpretations of non-GO MSigDB gene sets. The
corpus is deliberately diverse: C8 single-cell cell-type signatures, KEGG
legacy disease and pathway sets, Reactome and Hallmark sets, and a C7 cytokine
perturbation contrast. These exercise ontology-grounded `contexts` such as cell
type, tissue, disease, phenotype, experimental condition, and perturbagen.

The corpus also holds ranked signatures: a score for every measured gene,
such as rare-variant burden or a differential-expression contrast. See
[Ranked signatures](#ranked-signatures).

Each `genesets/<SET>.yaml` records, per GO term, a curator judgment
(role category, confidence, specificity), recovery status, enrichment stats,
and optional verbatim literature evidence. The data doubles as a
precision/recall gold standard for `genesets-rs` and as source content for a
future curated gene set browser.

## Layout
- `schema/genesets_interpretation.yaml` - LinkML schema (source of truth).
- `conf/` - OAK adapters and reference-validator config.
- `genesets/manifest.tsv` - index of curated sets (id, collection, biological context term, series, series_role).
- `genesets/<SET>.yaml` - one interpretation per gene set or ranked signature.
- `rankings/<NAME>.rnk.gz` - scores of ranked signatures (gzipped `gene`/`score` TSV).
- `cache/` - committed ontology label cache (written by linkml-term-validator).
- `references_cache/` - cited-paper cache (gitignored; rebuilt on demand).

## Workflow
1. `just curate-validate-schema` - sanity-check the schema's enum meanings.
2. Draft from enrichment output:
   `genesets-workflows curate draft MSIGDB:<SET> --enrichment-tsv <tsv> -o genesets/<SET>.yaml`
3. Adjudicate: set each association's `category`, `confidence`, `specificity`;
   add `curator_added` core terms the tool missed; add `evidence`.
4. `just curate-validate` - structural + term + reference + obsolescence +
   ranking validation.
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
insight-bearing term: **866 confirmatory, 87 mechanistic (~91% / 9%)** of 953.
The mechanistic minority concentrates in causal-gene phenotype sets (the
convergent molecular mechanisms behind a phenotype, e.g. RNA splicing / SMN-snRNP
assembly in motor-neuron disease), CRISPR-screen convergences (flavivirus host
factors → ER protein-biogenesis complexes; the *negative*-regulator brake modules
in T-cell and tumor-evasion screens), the down-arm of expression contrasts
(`JISON_SICKLE_CELL_DISEASE_DN`), and a few perturbations. Cell-type marker
signatures (the C8 single-cell sets) are almost entirely confirmatory — a marker
panel restates its own construction — so the insight value there is not the
mechanistic tag but the `membership_gap` calls where a fetal signature lacks its
cell type's mature effector machinery (e.g. catecholamine-synthesis enzymes
absent from fetal adrenal chromaffin cells; GAD/VGAT absent from fetal cortical
interneurons). This confirms the benchmark is predominantly an *annotation/method*
test bed; its insight-testing power is concentrated in a minority of sets.

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

## Ranked signatures

Some sources give a score for every measured gene rather than a list: a
differential-expression contrast, GWAS or exome gene-level statistics, a CRISPR
screen. Thresholding them throws information away, so they are curated as
ranked signatures, analyzed with `genesets-rs gsea`. What gets curated is the
same: per-GO-term judgments on the same axes. Four things differ.

- **The `ranking` block replaces membership.** It records the `source`, the
  exact `source_detail` and `download_url`, the `score_type`, whether the
  scores are `signed`, the `gene_id_type`, `n_ranked`, and the `file` and its
  `sha256`. The scores live in `rankings/<NAME>.rnk.gz`, a gzipped two-column
  TSV with a `gene`/`score` header that `genesets-rs gsea --ranks` reads as
  is. Each ranking is built from a recipe beside it, `rankings/<NAME>.recipe.yaml`:
  the source URL and its SHA-256, the table format, the gene and score columns,
  the id type, a transform (`identity`, `neg_log10`, `negate`), and a pinned
  HGNC release for mapping ids to approved symbols (previous and alias symbols
  are used only when unambiguous). `genesets-workflows curate build-ranking
  rankings/<NAME>.recipe.yaml` rebuilds the file byte for byte and reports how
  many source rows lacked a score or a symbol. Built files are committed only
  when the recipe sets `store: true`, because the source has no stable URL to
  rebuild from (the GWAS ATLAS MAGMA tables, served only through its web
  page); the rest are gitignored and rebuilt by `just curate-build-rankings`,
  which `just curate-validate` runs first. That also keeps publisher
  supplementary tables, often all-rights-reserved, out of the repository. The
  ranking is also the analysis universe: genes absent from it were not
  measured.
- **Direction is per association.** A signed ranking carries two programs, one
  at each pole, so each association records `direction` (UP toward high
  scores, DN toward low). An unsigned ranking, such as a -log10 p or a MAGMA Z,
  has only a high pole: its associations are all UP.
- **Ground judgments in the leading edge.** The leading edge (the genes that
  drive the score) is the ranked analogue of a set's membership. Read it, and
  name its genes in `curator_note`, as you would a set's genes.
- **`universe_gap` is the ranked `membership_gap`**: the genes that would carry
  a term were not measured. Genes that were measured but did not move are not a
  gap. That is a result, recorded in the description or a curator note.

`enrichment_stats` records which test produced it (`method`: `fisher`,
`gsea_ks` or `wilcoxon`; `algorithm`; `p_adjust_method`), with `es`, `nes`,
`set_size`, the leading edge and `peak_rank` for ranked tests. Choose the
statistic to suit the signal: GSEA's weighted KS for a concentrated top,
Wilcoxon for a modest shift across many genes, and say which in the
description. For scores derived from p-values (-log10 p, MAGMA), GSEA's
default `--weight 1` lets a few extreme genes swamp every running sum: on the
BMD MAGMA ranking, where -log10 p reaches 150, weighted KS finds no term and
unweighted KS (`--weight 0`) finds 181. Use Wilcoxon or `--weight 0` there.
For GWAS gene scores, also count independent loci in a leading edge: MAGMA
p-values of neighbouring genes are correlated through LD, so an HLA or APOE
block lights up whole terms as one finding.

CRISPR screens have their own traps, all visible in the curated screens:

- **Small complexes fall under the size floor.** The OST, signal peptidase
  and EMC complexes that define the flavivirus screens have 5-13 members, below
  `gsea`'s default `--min-size 15`; score them with `--min-size 5`.
- **Fitness.** Essential genes drop out of any survival or proliferation
  screen, and growth suppressors gain, whatever the selection. Record those
  terms as `nonspecific` (essential genes) or `false_association` (growth
  suppressors read as host factors).
- **Non-expression.** Genes the cells do not express escape the fitness cost,
  so whole-ranking tests shift olfactory receptors or heart genes. Record them
  as `false_association`.
- **Shared guides.** Clustered genes that share exons (PCDHA, UGT1A) are hit
  by the same guides and score as one: one finding, like an LD block.
- **Symbol hygiene.** Supplementary tables carry old symbols and Excel
  damage. `duplicates: prefer_exact` lets the row named by the approved symbol
  win; use `drop` when a collision may hide two genes (Excel-mangled MARCH/SEPT). Draft from a run with:

```bash
genesets-rs gsea --annotations gene_terms.tsv --terms terms.tsv --closure closure.tsv \
  --ranks rankings/<NAME>.rnk.gz --statistic wilcoxon --leading-edge-genes \
  --max-p-adjust 0.05 -o run.tsv
genesets-workflows curate draft-ranked LIT:<NAME> --gsea-tsv run.tsv \
  --ranks rankings/<NAME>.rnk.gz --source PMID:<id> --score-type <type> \
  --unsigned --method wilcoxon -o genesets/LIT_<NAME>.yaml
```

The draft fills the `ranking` block (count, path, checksum) and one
association per significant term, with `direction` from the sign of NES; for
an unsigned ranking it drops the meaningless low pole.

`curate validate` adds a ranking gate: the file exists and matches its
checksum and `n_ranked`, scores are well formed, an unsigned ranking has no
DN associations, every leading-edge gene is in the
ranking, `direction` agrees with the sign of NES, and ranked methods appear
only on ranked interpretations. `LIT:SCHIZOPHRENIA_SCHEMA_EXOME` is the first
example.

## Term obsolescence
`curate validate` adds a fourth gate: every ontology term id is swept against
its ontology's OAK `obsoletes()` set. This catches obsolete terms the
term-validator misses — the `sqlite:obo:*` builds retain labels for obsolete
classes, so an obsolete id paired with its old label still passes id+label
validation. An obsolete `id:` now fails the build (e.g. `GO:0050663` cytokine
secretion → use `GO:0032635`; `GO:0062023` → `GO:0031012`).
