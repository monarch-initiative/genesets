# Evidence-ablation eval against the curated benchmark

**Question:** How much of the curator-validated *core* GO biology does an
enrichment method recover when GOA human annotations are restricted by evidence
(all vs IBA-only vs IBA+IEA vs no-`contributes_to`)?

## Why this differs from `evals/go_iba_impact_*`

The `go_iba_impact` analyses compare annotation variants **to each other**
(retained / lost / gained calls) — there is no ground truth. This eval scores
each variant's enriched terms against **external ground truth**: the
curator-validated CORE terms in `curation/genesets/*.yaml`, bucketed by
`recovery_status`. The curated benchmark is what turns the impact analysis into
a real evaluation.

## Pipeline (in-repo tools)

1. **Membership** — fetch the benchmark's MSigDB sets' member symbols into a
   `queries.gmt`. mygeneset.info's `_id` is the MSigDB set name, so each curated
   `MSIGDB:<NAME>` set is fetched by querying `<NAME>` and taking the exact
   `_id`/`source=msigdb` hit. (11 MSigDB sets — mostly newer C8 single-cell
   clusters — are absent from mygeneset.) Most `LIT:` sets store markers in prose
   rather than full membership, but the 12 `LIT:GENETIC` sets are defined by
   short, explicit gene lists, so their membership was **captured directly from
   the primary sources**: GWAS / sequencing convergences (Satterstrom 2020 Table
   S2 [autism]; Trubetskoy 2022 prioritized + SCHEMA + C4 [schizophrenia]; the
   GWAS-Catalog loci for Jostins 2012 [IBD], Bellenguez 2022 [Alzheimer's], Nalls
   2019 [Parkinson's], van der Harst 2018 [CAD], Morris 2019 [bone density]) and
   CRISPR-screen hit tables (Daniloski + Wei [SARS-CoV-2]; Marceau/Zhang 2016
   flavivirus host-factor screens [OST/EMC convergence]; Shifrut 2018 T-cell
   proliferation regulators; Manguso 2017 in-vivo melanoma immunotherapy screen),
   plus one short curated mechanism panel (Bersuker/Doll 2019 FSP1-CoQ ferroptosis
   suppressors). All HGNC-normalized into `curation/genesets/lit_members.gmt` and
   folded into `queries.gmt`. Total evaluable: **133** (131 producing
   enrichment) — including 44 MSigDB C8 single-cell cell-type signatures
   (membership fetched from mygeneset.info). The two interferon sets added with
   the gene-set-informed factor-model docs are here too: the MSigDB pathway set
   `REACTOME_INTERFERON_SIGNALING` (fetched from mygeneset.info) and the six-gene
   `TYPE_I_INTERFERON_SCORE` (AGS-6 clinical panel, captured into
   `lit_members.gmt` — one of the few `LIT:DISEASE_ACTIVITY` sets with an explicit
   captured membership rather than prose markers). The five `LIT:PERTURBATION` RFX6-knockdown
   sets are likewise captured into `lit_members.gmt`: the Walker et al. (Nature
   2023) beta cell pair `RFX6_KNOCKDOWN_BETA_CELL_UP` / `_DN` (1,272 and 1,397
   members — the two largest queries in the eval, from Supplementary Data 2a),
   the 220-gene SCENIC regulon core of the same experiment
   (`RFX6_KNOCKDOWN_BETA_REGULON`, Supplementary Data 2c), and the Coykendall et
   al. (Diabetes 2024) alpha cell pair `RFX6_KNOCKDOWN_ALPHA_CELL_UP` / `_DN`
   (357 and 558 members, captured from the Supplementary Table S3/S4 PDFs; the
   parsed counts reproduce the paper's stated 357/558/915 exactly). Four of the
   five form `SERIES:RFX6_KNOCKDOWN_ISLET_CELL` with roles `beta_up`,
   `beta_down`, `alpha_up`, `alpha_down`.

2. **Annotation variants** —
   ```bash
   python3 scripts/prepare_go_eval.py --out-dir /tmp/iba_eval \
     --variants all,iba,iba_iea,no_contributes_to
   ```
   downloads GO + GOA, and writes per-variant `gene_terms.tsv` (evidence-filtered),
   shared `terms.tsv`/`closure.tsv`, a common all-GOA background, and a
   `genesets-rs matrix` config per variant (Bonferroni, `max_p_adjust 0.05`).
   Place the `queries.gmt` from step 1 at `/tmp/iba_eval/queries.gmt`.

3. **Enrichment** —
   ```bash
   for v in all iba iba_iea no_contributes_to; do
     genesets-rs run /tmp/iba_eval/$v/config.yaml
   done
   ```
   Each writes `<variant>/results.tsv` = significant (set, GO-term) hits.

4. **Score vs gold** —
   ```bash
   uv run --project python/genesets-workflows --extra curation \
     python scripts/score_method_vs_benchmark.py --eval-dir /tmp/iba_eval \
     --out /tmp/iba_eval/benchmark_scores.tsv
   ```

## Metrics

The PRIMARY metric scores the curated **biology** — recall of CORE terms
(`category` in `core_process`/`core_component`) — which is independent of the
gold's `recovery_status` labels. `recovery_status` is the curator's diagnostic
prediction and is used only descriptively; it is never refit to a method's
output (see "Guardrail").

- **`recall_core`** — recovered / total of all CORE terms. The headline: how
  much curator-validated biology a variant recovers.
- **`recall_supported`** — recall over the `annotation_supported` CORE subset
  (descriptive: of terms the curator *predicted* recoverable).
- **`gap_recovered`** — count of `annotation_gap` CORE terms a variant recovered.
  A **disagreement** between the curator's prediction and a tool run — a curator
  review item, not an automatic gold edit (see "Guardrail").
- **`unique_vs_baseline`** — supported-core terms a variant recovers that `all`
  does not.

## Headline result (2026, 133 evaluable / 131 scored sets, GOA `goa_human` current)

| variant | recall_core | gap_recovered (disagreements) |
|---|---|---|
| all | 0.559 (229/410) | 7 |
| no_contributes_to | 0.556 (228/410) | 7 |
| iba_iea | 0.437 (179/410) | 4 |
| iba | 0.334 (137/410) | 3 |

Refreshed as the five RFX6-knockdown sets were added. Each refresh re-scored the
same enrichment run against the previous gold to separate the new sets from GOA
drift, and both times the earlier figures reproduced exactly: 0.561 / 0.558 /
0.437 / 0.338 on 221, 220, 172, 133 of 394 core terms for the 128-set gold, and
0.565 / 0.562 / 0.440 / 0.341 on 227, 226, 177, 137 of 402 for the 130-set gold,
with identical gap_recovered counts throughout. The movement is therefore the
new core terms, not drift. Note the direction: recall_core *fell* from 0.565 to
0.559 because the RFX6 sets are hard — of their 16 core terms, all-GOA recovers
8.

1. **IBA carries ~2/3 of the core biology full GOA does** (recall_core 0.34 vs
   0.56); IEA recovers much of the difference (iba_iea 0.44).
2. **IBA is nearly a strict subset of all-GOA — it does not fill experimental
   gaps here.** Restricting to IBA loses 90 core terms and uniquely recovers
   only 2, both conserved-housekeeping cellular components (ribosome
   `GO:0005840`, nucleolus `GO:0005730`).
3. **The eval surfaces a review queue — it does not edit the gold.** 7
   `annotation_gap` core terms were recovered by standard all-GOA enrichment
   (e.g. SCHUHMACHER_MYC -> rRNA processing, TRAVAGLINI_CILIATED -> axoneme
   assembly, KEGG_RCC -> positive regulation of angiogenesis). These are
   *disagreements* between the curator's gap prediction and a tool run, queued
   for deliberate curator review against GOA facts — not auto-relabeled.

Calibration: even all-GOA recovers only ~58% of core terms under Bonferroni, so
a `core` + `annotation_supported` label means "the genes carry it", not "it
always reaches genome-wide significance".

## Confirmatory vs mechanistic (the insight split)

With the corpus-wide `insight` tags, recall splits by whether a term is
`confirmatory` (restates the set's construction) or `mechanistic` (a non-obvious
process — a genuine enrichment insight). Over the evaluable sets:

| variant | recall_confirm (n=659) | recall_mechan (n=86) |
|---|---|---|
| all | 0.563 | 0.349 |
| no_contributes_to | 0.561 | 0.349 |
| iba_iea | 0.431 | 0.244 |
| iba | 0.316 | 0.221 |

**Mechanistic insight is ~2x harder to recover than confirmatory biology** —
even all-GOA recovers only ~35% of mechanistic terms vs ~56% of confirmatory
ones. Standard enrichment surfaces the obvious and largely misses the
non-obvious convergent mechanisms the curators flagged (often the
`annotation_gap` ones).

The `LIT:GENETIC` sets (membership captured directly from the primary sources —
see Membership) make this concrete, per famous mechanism, under all-GOA. Among
the **GWAS / sequencing** convergences, two patterns emerge.

**EA recovers** the gene-dense, well-annotated convergences:
- Autism -> chromatin organization / transcription regulation (3/3)
- Bone mineral density -> Wnt signaling (2/2)
- SARS-CoV-2 -> V-ATPase / vacuolar acidification (3/6)

**EA is blind** to the celebrated mechanistic reframings — each invisible
because the convergence is spread thinly across many weakly-contributing genes
(and/or under-annotated):
- Alzheimer's -> microglia / lipid / complement / endocytosis: **0/5**
- Parkinson's -> lysosomal / autophagy: **0/3**
- Coronary artery disease -> vascular-wall ECM / NO signaling: **0/2**
- Schizophrenia -> complement-mediated synapse pruning (C4): **0/2**
- IBD -> autophagy: **missed**

The most field-reshaping insights of disease genetics — AD as a microglial/lipid
disease, PD as a lysosomal disorder, the schizophrenia C4 mechanism — are
exactly the ones standard enrichment cannot surface. A *curated* gold standard
captures them because a human asserts the mechanism from the literature, which
diffuse, many-gene, weak-per-gene enrichment never will. (The aggregate
`recall_mechan` ~0.34 hides this bimodality: the few recoverable mechanisms are
gene-dense and well-annotated; the celebrated ones are neither.)

### Why a mechanism is (in)visible — the CRISPR-screen sets

The four **CRISPR-screen** `LIT:GENETIC` sets sharpen the picture: unlike diffuse
GWAS loci, screen hits are direct functional convergences, yet they are recovered
very unevenly. Three distinct *reasons* a mechanism stays invisible emerge —
diffuseness (the GWAS story above) is only one of them:

- **Visible — the mechanism is a tight, well-annotated physical complex.**
  Flavivirus host-factor screens (Marceau/Zhang 2016) converge on ER protein-
  biogenesis machinery, and EA recovers **4/5** mechanistic terms: OST complex
  (`GO:0008250`), EMC complex (`GO:0072546`), ERAD (`GO:0036503`), ER
  (`GO:0005783`). The hits *are* the complex (STT3A/STT3B/RPN1/RPN2/OSTC →
  OST; EMC1-4/MMGT1 → EMC), so the convergence is gene-dense and densely
  annotated — the opposite of the GWAS case. (Only the 3-subunit signal
  peptidase complex is missed: too few member genes carry the annotation.)
- **Sign-invisible — the insight is the *direction* of regulation.** The two
  immune screens converge on *negative* regulators, and that sign is exactly
  what GO enrichment cannot see. Shifrut 2018 T-cell regulators: EA recovers the
  generic *T cell activation* / *TCR signaling* (confirmatory ✓✓) but misses
  **negative regulation of T cell activation** and **negative regulation of
  cytokine signaling** (0/2 mechanistic) — the brake module (CBLB, SOCS1,
  TNFAIP3, RASA2) is annotated to the activation processes, not to their
  repression. Manguso 2017 melanoma evasion: EA recovers the IFN-γ–sensing axis
  (*cellular response to type II interferon* ✓) but misses **negative regulation
  of IFN-γ signaling** — the PTPN2 insight, whose *loss* sensitizes tumors, the
  whole therapeutic point.
- **Size-invisible — the curated mechanism panel is too small to reach
  significance.** The 5-gene ferroptosis-suppressor panel (FSP1-CoQ axis,
  Bersuker/Doll 2019) recovers **0/5** terms, including the celebrated
  *ubiquinone biosynthetic process* mechanism: with n=5, no term clears
  Bonferroni regardless of annotation quality. A real, well-annotated mechanism
  can be invisible purely on set size.

So the gold standard now records three orthogonal failure modes for mechanistic
recall — **diffuse** (GWAS), **sign-blind** (directional regulation), and
**under-powered** (tiny panels) — only the first of which is about annotation
depth. The single recoverable family is the tight physical complex.

### A fourth failure mode — direction-of-change, from the RFX6 knockdown series

The five `LIT:PERTURBATION` RFX6 sets add a case the screens do not cover: the
arms of one perturbation disagree with the phenotype, and enrichment can only see
one of them. Four form `SERIES:RFX6_KNOCKDOWN_ISLET_CELL` (beta up/down from
Walker et al.; alpha up/down from Coykendall et al.), plus a regulon-level
granularity variant. Under all-GOA:

- **Beta DN — 5/5 core recovered.** Proteasome-mediated ubiquitin-dependent
  catabolism (`GO:0043161`), proteasome complex (`GO:0000502`), macroautophagy
  (`GO:0016236`), vesicle-mediated transport (`GO:0016192`) and cilium
  (`GO:0005929`) all clear Bonferroni — the tight-complex pattern again, with 25
  proteasome subunits, 10 ATG genes and 13 cilium/IFT genes literally in the set.
  The strongest of these (the ubiquitin-proteasome arm, LogP -39 in the paper's
  own Metascape output) is *not discussed anywhere in the paper*.
- **Beta UP — 1/3 core recovered.** The synaptic module (`GO:0099536`) is
  recovered; **insulin secretion** (`GO:0030073`) and **secretory granule**
  (`GO:0030141`) are missed, though ABCC8, SLC30A8, SLC2A2, G6PC2, CACNA1D,
  PCSK1/2, CPE, IAPP and MAFA are all members. Thirteen well-annotated genes in a
  1,272-gene query do not clear Bonferroni — and here the missed term is the
  whole point, since those transcripts go **up** while secretion is blunted.
- **Alpha DN — 1/4 core recovered, and the interesting miss is granularity, not
  power.** `GO:0030141` secretory granule is recovered at padj 9e-13 (55
  members), and the arm is thick with `regulation of secretion`, `regulation of
  insulin secretion` and `regulation of peptide hormone secretion` — yet the
  curated `peptide hormone secretion` (`GO:0030072`), `endocrine pancreas
  development` (`GO:0031018`) and `glucagon secretion` (`GO:0070091`) are all
  missed. The last was curated as an `annotation_gap` in advance and behaves
  exactly as predicted: SSTR2, RAPGEF4, KCNA5, CACNA1A, KCNMA1 and UCN3 are
  present but annotated only to generic hormone-secretion and ion-transport
  parents.
- **Alpha UP — 1/2 core recovered, and the same specificity effect.** `cellular
  response to stress` is recovered; `type I interferon-mediated signaling
  pathway` (`GO:0060337`) is not, even though MX1, IFIT1/2/3, ISG15, IRF7, STAT1
  and RIGI are members — the run instead returns the generic parents `innate
  immune response` (padj 4e-7) and `response to virus` (padj 5e-3). Specific,
  well-annotated terms lose to Bonferroni while their vague ancestors survive.
- **The precision items all fire.** `cell cycle` (`nonspecific` on beta DN,
  driven by proteasome subunits inside Reactome's APC/C steps), `brain
  development` (`marker_driven_plausible` on both beta UP and the regulon set)
  and `negative regulation of apoptotic process` (the corpus's first
  `false_association`, curated on alpha UP because that arm is a growth-arrest
  and pro-apoptotic-stress combination) are every one of them returned by
  all-GOA — four flagged-in-advance over-calls, four hits.

**A prediction this eval refuted.** The regulon set was added as a power control:
the 220-gene SCENIC core of the *same* beta cell experiment whose 1,272-gene UP
arm misses the secretory program. The authors' own analysis does recover it
there — KEGG Insulin secretion (LogP -7.9) and MODY (-9.2) via Metascape. The
expectation was that a smaller query would let GO-based enrichment find it too.
It does not: the regulon set returns 96 significant terms, all of them neuronal
(synapse padj 1e-18, presynapse, synapse assembly) and **not one secretory term
at all**. The gap is not query size but ontology — NKX6-1, SLC30A8, INS, PCLO,
CACNA1C/D, KCNB2 and KCNMB2 sit inside KEGG's channel-centric insulin-secretion
pathway while being thinly annotated to GO's insulin secretion. A GO-only
benchmark cannot see a pathway that KEGG draws and GO does not.

### One program at three membership granularities — the interferon gradient

The corpus curates ONE biological program — the type I interferon response — at
three deliberately different membership sizes (the 6-gene clinical score, the
~20-gene discovery signature, and the 296-gene Reactome pathway; see the
gene-set-informed factor-model section of the Competitive Landscape). Two of the
three carry captured membership in `queries.gmt`, and they sit at the extremes,
so the eval reads out the size effect cleanly:

- **`REACTOME_INTERFERON_SIGNALING`** (296 genes) recovers **4/4** of its core
  terms under all-GOA — type I and type II interferon-mediated signaling, antigen
  processing and presentation, and the MHC protein complex — because the pathway
  set is gene-dense and exhaustively annotated across every sub-program.
- **`TYPE_I_INTERFERON_SCORE`** (the 6-gene AGS-6 clinical panel) recovers **0/3**
  core terms. Two are curator-labelled `membership_gap` by construction (antigen
  presentation and the JAK-STAT signaling machinery — carrier genes deliberately
  absent from a minimal ISG-effector readout), so they are correctly unrecovered;
  the third, *cellular response to type I interferon*, is `annotation_supported`
  (all six genes carry it) yet still misses Bonferroni significance at n=6. It is
  the interferon twin of the ferroptosis panel: a real, well-annotated program
  that is **size-invisible** to enrichment.

Same biology, opposite ends of the completeness axis — a compact demonstration
that `annotation_supported` means "the genes carry it", not "it clears
significance", and a curated fixture for the membership-refinement question these
factor models address from data.

### Single-cell C8 signatures — when a cell type's defining function is absent
from its own signature

The 36 MSigDB C8 single-cell cell-type signatures add a different lesson. Several
of the fetal **Descartes** signatures are dominated by lineage-*specification*
transcription factors and lincRNAs, with the cell type's textbook **mature
effector machinery absent** — so its defining function is a curator-asserted
`membership_gap`, and the eval confirms it is unrecovered, exactly as predicted:

- Fetal cerebrum inhibitory neurons → **GABA biosynthesis / GABAergic
  transmission**: `membership_gap` (GAD1/GAD2/SLC32A1 absent; the set carries the
  DLX/ARX/LHX8 specification program instead). Unrecovered. ✓predicted
- Fetal adrenal chromaffin cells → **catecholamine biosynthesis** and the
  **chromaffin granule**: `membership_gap` (TH/DBH/DDC/PNMT and CHGA/CHGB absent;
  the set carries the sympathoadrenal program SLC6A2/ISL1/HAND1/NTRK1).
  Unrecovered. ✓predicted

This is the gold doing its job in the opposite direction from the GWAS sets: not
"a real mechanism the tool misses", but "the tool *correctly* fails to find a
function whose genes aren't present" — a curator can only flag that by knowing
what the signature *should* contain but doesn't.

Among the adult/crisp signatures, recovery splits the usual way: effector-dense
sets recover ~half their core (microglia 3/6 — phagocytosis, activation,
inflammatory response; intestinal epithelium 3/6 — brush border, digestion,
xenobiotic; photoreceptor 3/6 — the structural outer-segment/cilium terms, while
the fetally-incomplete phototransduction cascade is not). Two crisp sets recover
**0** core despite every term being `annotation_supported` (pancreatic delta 0/6;
airway goblet 0/5): the carrier genes are a minority of a large signature and
never clear Bonferroni — the calibration point that `annotation_supported` means
"the genes carry it", not "it reaches genome-wide significance".

A second C8 batch (smooth muscle, basophil/mast, megakaryocyte, RPE,
enteroendocrine, Schwann, chondrocyte, placental trophoblast) widens the
recoverability spread further. The two extremes:
- **Megakaryocyte 6/6** — the most completely recovered set in the corpus:
  platelet formation, aggregation, alpha granule, MK differentiation, coagulation,
  hemostasis. Platelet biology is gene-dense *and* exhaustively GOA-annotated, so
  every core term clears significance. Chondrocyte (5/6) and Schwann (4/6, the
  PNS-specific myelination term an `annotation_gap`) follow.
- **RPE 1/6, mast 1/6, enteroendocrine 1/5** — `annotation_supported` cores that
  largely miss significance (small or dilute signatures), plus more correctly-
  unrecovered gaps: RPE melanin biosynthesis (`membership_gap` — TYR/TYRP1/MITF
  absent from fetal, not-yet-pigmented RPE) and trophoblast placenta-development /
  syncytin fusion (`membership_gap` — the fusogens ERVW-1/ERVFRD-1/GCM1 absent).

These five `SERIES:` contrast groups resolve cleanly — alveolar_type_1↔type_2,
islet alpha/beta/delta, cerebrum excitatory↔inhibitory, myelinating glia
central↔peripheral (oligodendrocyte vs Schwann), and muscle smooth/skeletal/cardiac
— so the eval can check that the poles map to contrasting GO interpretations.

## Guardrail: the eval must not refit the gold

The eval is a measurement, not an editor. The gold's `category` (biology) is the
ground truth scored here and stays the authority. `recovery_status` is the
curator's diagnostic prediction; it is checked against GOA facts during
curation, never auto-updated to match whatever a given method + annotation
snapshot + p-threshold happens to recover. Refitting `recovery_status` to a tool
would make recall-vs-gold circular — the tool could not be "wrong" because we'd
have redefined truth to match it — and would forfeit the gold's independence
across methods and GOA versions. The `gap_recovered` disagreements are review
*inputs*, adjudicated deliberately on the merits, not automatic relabels.

Generated tables (`/tmp/iba_eval/`) are not committed; rerun the pipeline to
regenerate. The scorer is `scripts/score_method_vs_benchmark.py`.
