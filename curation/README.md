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
- `ontology_migrations.tsv` - audit log of every ontology-drift migration (see "Ontology drift").
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
4. `just curate-validate` - structural + term + reference + ontology-drift +
   ranking validation (reports every failing file). After a GO release, run
   `just curate-drift` / `just curate-migrate` (see "Ontology drift").
5. `just curate-report` - precision/recall/F1.

## Categories (biology)
`core_process`, `core_component`, `supporting_process`,
`marker_driven_plausible`, `nonspecific`, `false_association`.

The `category` is an assertion about the **biology** and is *not* driven by the
current state of GO annotations. A term that is biologically core stays core
even when no annotation supports it — we never downgrade truth to match
incomplete data.

## Noise labels and precision

`nonspecific` (true but uninformative for this set) and `false_association`
(an artefact, or wrong for the set) are the noise categories. Give each noise
term an `artifact_cause`, which says why it is noise:

| `artifact_cause` | meaning |
|---|---|
| `multifunctional_gene_bias` | driven by a few highly annotated genes; disappears under a multifunctionality correction (Ballouz et al. 2017, PMID:28204549) |
| `over_annotation` | carried by annotations that are themselves wrong or over-propagated |
| `housekeeping_contamination` | driven by housekeeping or ribosomal genes that entered the set for technical reasons |
| `overly_general_term` | true of the genes but too general to say anything about the set |
| `screen_or_assay_artifact` | fitness genes, shared guides, batch effects and the like |
| `wrong_biology` | wrong for the set, with none of the causes above |

**Scoring precision needs every hit labelled.** A curated set normally lists
a handful of terms, so a reported hit with no association is ambiguous: it is
unreviewed, not necessarily noise. To make precision scoreable for a set:

1. Run every annotation variant you want to compare (e.g. `all`, `iba`,
   `iba_iea`, `aigr_prune`, `aigr_full`) and take the **union** of their
   significant hits for the set (pooling). A variant that surfaces new noise
   is then penalised for it.
2. Label every pooled hit, **blind to which variant reported it**. Judge the
   term, its size and the member genes that carry it; never "this variant is
   expected to be better".
3. Record the pool in the set's `hit_adjudication` block (variants, test,
   cutoff, correction scope, GOA/GO release, date).

The validator then requires a category on every association of that set, and
the eval reports per-variant precision (`signal / (signal + noise)`, with
noise broken down by `artifact_cause`) over adjudicated sets only. Terms the
source authors call artefacts can also be added as curator-added
`false_association` terms, so that any variant reporting them is penalised.
The Ballouz 2017 case studies (`LIT_GILMAN_2011_AUTISM_CNV`,
`LIT_MANALO_2005_HYPOXIA_HIF1_UP`, `LIT_SCHMIDT_KASTNER_2012_SCHIZOPHRENIA`;
sources in `curation/sources/ballouz2017/`) are the first adjudicated sets.

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

A third literature collection, `LIT:PERTURBATION`, is the `LIT:` analogue of
MSigDB's `C2:CGP`: the transcriptional **readout of a deliberate perturbation**
of one gene, captured from the defining paper. It is neither driver nor disease
activity — the perturbed gene is the *context* (a `genetic_perturbation`
`BiologicalContext`, typically a `PR:` protein term), and the members are the
downstream response. Recovery-status defaults follow the activity pattern
(`annotation_supported`), since the differentially-expressed members carry the
annotations for the processes they read out. The RFX6-knockdown pair
(`LIT:RFX6_KNOCKDOWN_BETA_CELL_UP` / `_DN`, Walker et al. Nature 2023) is the
worked example, and it shows why the perturbed gene must not be confused with
the membership: RFX6 itself is not a member of either arm.

A fourth, `LIT:PHENOTYPIC_PROFILING`, holds **connectivity (guilt-by-similarity)
sets** from image-based or other phenotypic profiling: the members are the
annotated targets of reference compounds whose profiles match a query
perturbagen, so they are the genes whose inhibition *phenocopies* the query.
The query's own target is the `chemical_perturbagen` context, not a member, and
the pathway it drives is typically a `membership_gap`. These sets carry
construction biases worth recording as negative controls. Annotated compound
libraries are dominated by kinase inhibitors, so kinase/phosphorylation terms
are `nonspecific`. A morphological match can also come from a shared cell shape
rather than a shared mechanism, which makes such co-matches at most
`marker_driven_plausible`. The worked example is
`LIT:PPARG_INVERSE_AGONIST_PHENOCOPY_TARGETS`: the JUMP Cell Painting compound
targets that phenocopy covalent PPARG inverse agonists (CDK4/6, RRM1/2, TOP2A,
…), which resolve to the CDK4/6-Rb G1/S checkpoint and senescence.

### Model-organism sets
A set is curated in the species it was defined in (`taxon`), and must be scored
against that species' GAF, never `goa_human`. Three non-human sources are in use:
native-mouse MSigDB M2 perturbation sets (`MSIGDB:` ids, `collection: M2:CGP`,
`msigdb_release: v2025.1.Mm`), and CTD chemical-gene interaction sets for worm,
fly and zebrafish as served by mygeneset.info (`MYGENESET:<MeSH>_<taxid>`,
`collection: CTD:chem_gene_ixns`; context is the chemical's `CHEBI:` or `MESH:`
term); and yeast genome-wide deletion screens as `LIT:GENETIC` sets, whose
membership is SGD's curated null-mutant phenotypes (budding yeast,
`taxon: NCBITaxon:559292`) or PomBase's curated deletion FYPO phenotypes (fission
yeast, `NCBITaxon:4896`) from the one defining paper. CTD sets pool mixed-direction literature (sensitivity phenotypes,
expression, proteomics), so they are `direction: NA` and have no single defining
paper to quote. Membership for every non-human set is snapshotted in
`genesets/model_organism_members.gmt`. Ground and adjudicate `recovery_status`
on the species' GOA (MGI, WormBase, FlyBase, ZFIN, SGD), which often has different
gaps from human GOA. For example, mouse Saa1, Saa2 and Hp are not annotated to
acute-phase response. Eval: `evals/model_organism_benchmark/`.

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

## Published analyses (prior enrichment results)

`evidence` records the mechanisms a paper asserts. `published_analyses`
records what enrichment analyses of the set returned in the literature: the
defining paper's own run (`analysis_role: original`) and later reanalyses by
others (`reanalysis`). Each entry gives the reference, tool, method,
`annotation_source` and `annotation_date`, the `membership_size` it ran on
with a `membership_note` relating it to this entry's membership, the
`correction` and `significance_threshold`, `n_significant` with its
`significant_unit` (GO terms, or PANTHER's groups of related terms), a
`data_url` for deposited inputs or results, a `summary`, and an optional
verbatim `snippet`, which the reference validator checks like an evidence
snippet.

`reported_terms` lists the terms the paper's argument rests on, with the
stats as printed (`p_value`, `p_adjust`, and `p_adjust_corrected` when the
paper also reports a bias-corrected value, e.g. Ballouz et al.'s
multifunctionality correction), whether each was `significant`, and the
`authors_stance` (`highlighted` / `reported` / `discounted`). Papers from
older GO releases cite ids that have since been renamed or retired, so
`reported_id` and `reported_label` are copied as printed and not validated;
`term` is the curator's mapping to a current GO term, which goes through the
id+label gate and the ontology-drift gate (and is rewritten by
`just curate-migrate` like any other term; the `reported_*` lines are never
touched). Leave `term` out when nothing current fits.

Keep one entry per paper and setting: the same list analysed with all GO
annotations and with PAN-GO alone is two entries. Record membership
differences explicitly. Reanalyses often run on a different version of the
list (an older MSigDB mapping, a larger cluster from the same study), and
results are not comparable without that context.

Published analyses are a historical record, not gold. They are never scored
and never set `category`, `recovery_status` or `insight`. A curator may cite
them in a `curator_note`, and should check whether a reported result still
holds under current GOA. For example, the peptidyl-proline modification
result in `LIT_MANALO_2005_HYPOXIA_HIF1_UP` has no carriers today, because the
prolyl hydroxylases are now annotated with molecular-function terms instead.
The three case studies of Ballouz et al. (NAR 2017) and the PAN-GO
functionome paper (Nature 2025) carry this block: `LIT_GILMAN_2011_AUTISM_CNV`,
`LIT_MANALO_2005_HYPOXIA_HIF1_UP` and `LIT_SCHMIDT_KASTNER_2012_SCHIZOPHRENIA`
(sources in `curation/sources/ballouz2017/`).

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

## Ontology drift (obsolescence, merges, relabels)
GO and the other cited ontologies change underneath the gold standard. Each
change costs something different. An obsolete term can never be recovered, because GOA stops
annotating to it, so a core term that went obsolete silently lowers recall. A
merged id stops resolving. A relabel passes the term-validator, which accepts
synonyms, but leaves the file out of date. The repo treats these changes as routine
maintenance with an audit trail, not as surprise build failures:

- **Gate.** `curate validate`'s drift check (gate 4) classifies every context and
  association term against the current OAK ontology. Obsolete, merged and
  unknown ids are errors; a relabel of a live term is a warning, so a GO rename
  never breaks the build. Each message names its fix.
- **Report.** `just curate-drift` scans the whole corpus, prints each ontology's
  version, and exits 1 while anything is unresolved (usable as a CI check).
- **Migrate.** `just curate-migrate` rewrites what the ontology itself declares:
  relabels, merges (alt ids → primary) and obsoletions with exactly one live
  `replaced_by`. Only the term's `id:`/`label:` lines change, and each rewrite is
  re-read and verified. A fix whose target is already asserted in the same file
  is held back as `needs_review` so associations are never silently duplicated.
- **Curator decisions.** An obsoletion with only `consider` suggestions or no
  replacement needs a curator:
  `genesets-workflows curate drift FILE --map OLD=NEW --note "why"`. Choose
  from the members' actual annotations, then update the `curator_note`.
- **Audit log.** Every applied change is appended to
  `curation/ontology_migrations.tsv` (old/new id and label, kind, ontology
  version, GO tracker issue, note).

Migration never edits `category`, `recovery_status` or `insight`. It follows
the ontology's declarations and keeps the curator's judgment. When the
replacement is narrower or broader, re-check `recovery_status` against GOA as
part of the review. For example, on 2026-11-14 no CAD member was still annotated
to obsolete `GO:0007263`, while four carry its replacement `GO:0038060`.

The eval side pins and records its GO version: `prepare_go_eval.py
--go-release YYYY-MM-DD` uses a dated release, resolving that release's archive
layout, and every run records the OBO `data-version` in `metadata.json`.
`score_method_vs_benchmark.py` prints that version and warns, or fails with
`--strict`, when a gold GO term is obsolete or unknown in the eval's release,
so version skew between the gold and the eval is never mistaken for a recall
change.

## Domain & Multi-Species Sub-Corpora
- [`curation/cgsa/`](cgsa/README.md) — Curated GO interpretations for the Hu et al. 2024 Nature Methods LLM gene set evaluation panel (`evals/cgsa/`).
- [`curation/plants/`](plants/README.md) — Multi-species plant GO-interpretation gold standard (30 gene sets across 6 vascular plant species) spanning *Arabidopsis thaliana* (`NCBITaxon:3702`) and five DOE/BER bioenergy & oilseed species (*Populus trichocarpa* `NCBITaxon:3694`, *Panicum virgatum* `NCBITaxon:38727`, *Sorghum bicolor* `NCBITaxon:4558`, *Thlaspi arvense* `NCBITaxon:13288`, and *Brachypodium distachyon* `NCBITaxon:15368`), paired with the multi-species GOA evidence variant benchmark in [`evals/plants/`](../evals/plants/README.md).
