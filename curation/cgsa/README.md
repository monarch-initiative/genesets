# cGSA-derived curated interpretations

Curated GO interpretations of gene sets drawn from the **cGSA benchmark** (Wang
et al., *Knowledge-guided contextual gene set analysis with large language
models*, Bioinformatics 2026, btag214; https://github.com/ncbi-nlp/cGSA). These
are kept in a **separate corpus** from `curation/genesets/` on purpose — see
"Why isolated" below.

The source data is ingested under `evals/cgsa/` (see that directory's README).
Each set here is one differentially-expressed-gene (DEG) module from the
benchmark, curated with the same discipline as the native corpus
(`curation/README.md` is the authority): every association grounded in the
**actual** DEG membership, every ontology id OLS-verified and swept for
obsolescence, biology-first `category`, orthogonal `recovery_status`, and
`insight`.

## Layout
- `genesets/CGSA_<SET>.yaml` — one interpretation per DEG module.
- `genesets/cgsa_members.gmt` — the module memberships (from the benchmark).
- `genesets/manifest.tsv` — index (same columns as the native manifest).

## Current sets

Thirty-two DEG modules across twenty studies, spanning melanoma, glioblastoma,
diabetic cardiomyopathy, Tfh-cell regulation, NASH liver, pancreatic cancer,
osteoporosis, squamous cell carcinoma, rheumatoid arthritis, brain,
Mycobacterium avium infection, leukemia, Alzheimer disease, nasopharyngeal
carcinoma, Williams-Beuren syndrome, Parkinson disease, tumor-associated
macrophages, intervertebral-disc (nucleus pulposus), and neural-stem-cell
development, in both human and mouse — covering antigen presentation, synaptic
organization, chromatin/nucleosome, muscle contraction, lipid/steroid metabolism,
Rho-GTPase actin remodeling, nucleolar ribosome biogenesis, elastin/basement-
membrane matrix, ciliary axoneme, centrosome/microtubule, and mitotic
chromosome segregation. Four studies are represented as multi-module series
(melanoma, glioblastoma, NASH, diabetic cardiomyopathy), and several studies now
have complementary modules curated (Parkinson synapse + cilia; pancreatic Golgi +
Rho-actin + chromosome).

Instructive "what the membership actually supports" calls recur: the NASH
ribosome module records KEGG neurodegeneration pathways as `false_association`
artifacts; the glioblastoma-stress module downgrades a broad "response to stress"
ground truth to `nonspecific`; the Williams-Beuren module curates the
elastin/ECM signal (disease-relevant) and marks a strong olfactory-receptor
enrichment `nonspecific` as a gene-family artifact; and the rheumatoid-arthritis
"T cell" set is curated to the striated-muscle biology its membership carries
(documented mismatch, see issue #9).

> Most of these were drafted with **`evals/cgsa/draft/draft_cgsa.py`**, which
> does the deterministic gene-list work — real genesets-rs enrichment against
> species-matched GOA, so each candidate GO term arrives with its true carrier
> genes and stats, plus a GT/membership mismatch flag — and leaves the GO-term
> `category`/`insight` and context grounding to the curator. See
> `evals/cgsa/draft/README.md`. All 71 remaining benchmark sets have been drafted
> and triaged in **`evals/cgsa/draft/triage.tsv`** (53 curatable, ranked by
> informative-term richness; empty/thin sets flagged for skip). The drafter also
> caught a benchmark mislabel: PMC9205785-2, ground-truthed as a "T cell receptor
> complex" set, is actually a striated-**muscle** signature (sarcomere/myofibril)
> — curated to the biology its membership supports, with the mismatch documented.

**PMC10202813 melanoma modules** — one study whose DEGs the benchmark splits into
modules, linked as `SERIES:CGSA_PMC10202813_MODULES` to show that distinct modules
of one study resolve to **contrasting** GO interpretations:

| set | series_role | defining biology | grounding note |
|-----|-------------|------------------|----------------|
| `CGSA_PMC10202813_3_INTERFERON` | interferon_response | type I/II IFN response, antiviral defense | clean ISG panel (OAS/MX/IFIT/ISG15) — `annotation_supported`, confirmatory |
| `CGSA_PMC10202813_6_PIGMENT` | pigmentation_and_proliferation | melanin biosynthesis + melanosome | pigment core is **thin** (OCA2), co-mixed with a proliferation signal (CDC45/AURKB/PBK) |
| `CGSA_PMC10202813_2_SENESCENCE` | senescence_and_immune | cellular senescence + T-cell activation | two coherent arms (CDKN2A/TP53; PTPRC/TIGIT/SELL) plus melanoma drivers |

**PMC9232499 glioblastoma modules** — paired arms of one study, linked as
`SERIES:CGSA_PMC9232499_MODULES`:

| set | series_role | defining biology | grounding note |
|-----|-------------|------------------|----------------|
| `CGSA_PMC9232499_1_GBM_MIGRATION` | migration_ecm | cell migration + ECM remodeling | clean (Mmp9/Mmp10/Col1a1/Adamts2/Cyr61) |
| `CGSA_PMC9232499_2_GBM_STRESS` | stress_defense | innate/IFN defense response | GT's broad "response to stress" recorded `nonspecific`; the informative arm is Gbp1/Ciita/Fgl2 |

**Standalone modules:**

| set | defining biology | grounding note |
|-----|------------------|----------------|
| `CGSA_PMID30017245_MELANOMA_MITOSIS` | mitotic cell cycle | 12/12 canonical mitotic markers (CDK1/CCNB1/AURKB/PLK1/…) |
| `CGSA_PMC8831505_1_TFH_HYPOXIA` | hypoxia + glycolysis | HIF/glycolytic block (Slc2a1/Ldha/Hk2/Pdk1/Bnip3) in Tfh cells |
| `CGSA_PMC10761883_2_DIABETIC_CM_HIF` | hypoxia/HIF + PPAR/lipid | HIF targets (Bnip3/Adm/Angptl4/Vegfa) + glycolysis; KEGG "HIF-1 signaling" → GO response to hypoxia |
| `CGSA_PMC10156590_2_PANC_GOLGI` | ER–Golgi vesicle transport | COPI coatomer (Copb1/Copb2/Copg1) + Arf1/Rab1a; pancreatic cancer |
| `CGSA_PMC8983726_1_OSTEOPOROSIS_BONE` | ossification / osteoblast differentiation | Alpl/Bmp2/Bmp4/Msx2/Satb2; osteoporosis (drafted semi-automatically) |
| `CGSA_PMC9750880_1_SCC_MIGRATION` | actin/migration + angiogenesis | Abl2/Arpc1b/Bcar1/Fermt2 + Col4a1/Ackr3; squamous cell carcinoma (drafted semi-automatically) |

**PMC10772820 NASH-liver modules** — one mouse steatohepatitis study split into
four KEGG modules, linked as `SERIES:CGSA_PMC10772820_MODULES`:

| set | series_role | defining biology | grounding note |
|-----|-------------|------------------|----------------|
| `CGSA_PMC10772820_1_NASH_COMPLEMENT` | complement_coagulation | complement + coagulation cascade | hepatic C3/Cfb + F2/Fga/Fgb/Fgg/Plg |
| `CGSA_PMC10772820_3_NASH_RIBOSOME` | ribosome_translation | cytoplasmic translation | 12/12 Rpl/Rps; KEGG "Parkinson/Huntington/OXPHOS" ground truth recorded as `false_association` artifacts |
| `CGSA_PMC10772820_2_NASH_ECM` | ecm_fibrosis | ECM organization / fibrosis | Col1a1/Col4a1/Col6a1 + Lama2/Lamb1 |
| `CGSA_PMC10772820_5_NASH_CHEMOKINE` | chemokine_inflammation | chemokine signaling / inflammation | Ccl2/Ccl4/Ccl5/Cxcl9/Ccr2 |

The NASH ribosome module is a second instructive case (alongside GBM_STRESS): its
ground truth includes KEGG neurodegeneration/OXPHOS disease pathways that merely
re-list ribosomal and mitochondrial genes — recorded as `false_association`, not
curated as disease biology present in a liver DEG set.

Several modules were chosen precisely because they are **mixed** DEG clusters:
they exercise the membership-grounding rule (the paper's top ground-truth term
can be carried by very few set members), which is the within-cluster heterogeneity
cGSA's community detection is built to handle.

### Curation caveat: ground-truth vs membership mismatch

Not every benchmark set is curatable. Several sets' paper-derived ground-truth
pathways do **not** match their released DEG membership — e.g. `PMC11285963-1`
(pancreatic cancer) carries an epithelial/tumor membership (KRT/CEACAM/MUC1) but
its ground truth is entirely immune (T-cell activation, neutrophil
degranulation), and the two rheumatoid-arthritis lymphocyte sets list B-/T-cell
ground truth while only 1–2 canonical lymphocyte markers are actually present.
Curating those against the ground truth would violate the ground-in-membership
rule, so they are **skipped**. Every set here was screened for marker support
before curation.

## Why isolated from `curation/genesets/`

The native gold is scored by `scripts/score_method_vs_benchmark.py`, which globs
`curation/genesets/*.yaml` and measures recall against **human GOA** enrichment.
The cGSA benchmark is 52%% mouse and consists of large (~800-gene) DEG contrasts;
dropping those into the human-GOA precision/recall gold would silently distort
`recall_core`. Keeping them here preserves the native benchmark's integrity
while still subjecting every cGSA-derived YAML to the same 4-gate validator.

## Validate

```bash
# same 4 gates as the native corpus (structural, term id+label, reference, obsolescence)
uv run --project python/genesets-workflows --extra curation \
  genesets-workflows curate validate curation/cgsa/genesets/CGSA_PMC10202813_3_INTERFERON.yaml
```

> Environment note: the validator's OAK/eutils import chain needs
> `pkg_resources`, which setuptools ≥ 81 removed. If validation dies with
> `ModuleNotFoundError: No module named 'pkg_resources'`, pin `setuptools<81` in
> the workflows venv (already declared in the `curation` extra).
