# Does applying AI Gene Review to GOA improve enrichment of external gene sets?

**Question:** if the review decisions in
[ai-gene-review](https://github.com/ai4curation/ai-gene-review) (AIGR) were
applied to GOA human (over-annotations removed, terms corrected, new terms
added), would GO enrichment of *external* gene sets get better?

The test sets are external gene sets, not GO-derived ones:

| arm | sets | gold | what "better" means |
|---|---|---|---|
| **curated** | the 128 evaluable benchmark sets (MSigDB H/C2/C8 + `LIT:` GWAS/CRISPR) | CORE terms in `curation/genesets/*.yaml` | recover as many CORE terms, with the same or stronger p-values |
| **chem** | 270 MSigDB C2:CGP chemical-exposure signatures (`description:PubChem`, 10–1000 genes) | none | fewer hits, and the hits removed are uninformative |
| **onco** | 189 MSigDB C6 oncogenic signatures | none | as for chem |

The benchmark sets are also run through the no-gold metrics as `benchmark`.

## Variants

Built by `scripts/prepare_aigr_changeset_eval.py` from the AIGR change set
(`scripts/export_goa_changeset.py` in ai-gene-review). Edits are applied **row
by row** to the current `goa_human.gaf`, matching on (UniProt accession, GO
term, evidence code, reference). A (gene, term) pair survives if any of its
supporting rows survives.

| variant | edit |
|---|---|
| `all` | current GOA, unchanged (baseline; also defines the common background) |
| `aigr_remove` | drop REMOVE rows |
| `aigr_prune` | drop REMOVE + MARK_AS_OVER_ANNOTATED rows — **the main arm** |
| `aigr_full` | `aigr_prune` + MODIFY → replacement term(s) + NEW annotations |
| `aigr_core` | reviewed rows kept only if ACCEPT (+ MODIFY replacements, NEW): also drops KEEP_AS_NON_CORE and UNDECIDED |
| `null_prune_00..19` | negative control: from the same reviewed rows, drop a random sample with `aigr_prune`'s per-evidence-code row counts |

The null is what makes the result interpretable: any removal of annotations
shrinks enrichment output, so the question is whether **AIGR's** removals do
something different from removing the same number of rows at random.

## Pipeline

```bash
E=/tmp/aigr_eval
# 1. change set (in ai-gene-review)
uv run python scripts/export_goa_changeset.py --organism human -o exports/goa_changeset_human.tsv
# 2. variants + nulls (downloads GO + GOA)
python3 scripts/prepare_aigr_changeset_eval.py --out-dir $E \
  --changeset ../ai-gene-review/exports/goa_changeset_human.tsv
# 3. query sets
python3 scripts/fetch_benchmark_queries.py --out $E/queries.gmt
cd python/genesets-workflows
uv run python -m genesets_workflows.sources.mygeneset --query 'msigdb.subcategory.code:CGP AND description:PubChem' \
  --out-dir $E/breadth/chem --source-filter msigdb --min-genes 10 --max-genes 1000
uv run python -m genesets_workflows.sources.mygeneset --query 'msigdb.category.code:C6' \
  --out-dir $E/breadth/onco --source-filter msigdb --min-genes 10 --max-genes 1000
mkdir -p $E/breadth/benchmark && cp $E/queries.gmt $E/breadth/benchmark/
# optional: the rest of human MSigDB except C5 (GO-derived), ~16k sets
f(){ uv run python -m genesets_workflows.sources.mygeneset --query "$2" --out-dir $E/breadth/$1 \
       --source-filter msigdb --min-genes 10 --max-genes 1000; }
f hallmark 'msigdb.category.code:H'
f cgp_other 'msigdb.subcategory.code:CGP AND NOT description:PubChem'
f cp_pathways 'msigdb.subcategory.code:("cp:biocarta" OR "cp:pid" OR "cp:wikipathways" OR "cp:kegg_legacy" OR "cp:kegg_medicus")'
f reactome 'msigdb.subcategory.code:"cp:reactome"'
f c3_regulatory 'msigdb.category.code:C3'
f c4_cancer_modules 'msigdb.category.code:C4'
f c7_immunologic 'msigdb.category.code:C7'
f c8_celltype 'msigdb.category.code:C8'
# 4. run genesets-rs for every variant x arm, and score
cargo build --release
uv run --extra curation python ../../scripts/run_aigr_changeset_eval.py --eval-dir $E
```

Outputs are in `$E/report/`:
- `curated_summary.tsv`
- `curated_per_set_deltas.tsv`
- `curated_core_term_changes.tsv`: each CORE term gained or lost, with the set members whose annotation changed.
- `breadth_summary.tsv`
- `report.json`: includes the paired tests and the null comparison.

Generated files are not committed.

## Result (2026-10-03; AIGR `4b945bf`, 2,215 human reviews; GOA `goa_human` current)

All numbers below use the default inputs: direct `protein binding`
(GO:0005515) annotations are dropped from every variant, baseline included.
The section "With protein binding kept" compares this with the earlier run.

**Change-set coverage.** 105,086 reviewed (non-NOT) annotation rows; 98,195
(93%) match a row in current GOA, 82,990 of them outside `protein binding`.
1,470 NEW rows. Reviewed genes make up a median **17%** of the members of a
benchmark set (16% chem, 12% onco).

Row removals mostly do not remove pairs. REMOVE + MARK_AS_OVER_ANNOTATED
drop 5,417 non-`protein binding` rows, but only 4,324 (gene, term) pairs.
In most cases another row, unreviewed or accepted, still asserts the same
pair, or a kept descendant term still implies it.

### 1. Recall of the curated biology: no gain, a small loss at the edge of the null

| variant | recall_core (125 sets) | recall_confirm | recall_mechan | significant terms |
|---|---|---|---|---|
| all | 0.565 (221/391) | 0.568 | 0.351 | 28,915 |
| aigr_remove | 0.563 (220) | 0.568 | 0.351 | 28,841 |
| **aigr_prune** | **0.558 (218)** | 0.565 | 0.351 | **28,062** |
| aigr_full | 0.560 (219) | 0.562 | 0.338 | 28,001 |
| aigr_core | 0.535 (209) | 0.531 | 0.312 | 23,595 |
| null_prune (mean of 20) | 0.562 (219.6; range 218–221) | 0.565 | 0.348 | 28,550 (28,483–28,615) |

- `aigr_prune` loses 3 CORE terms:
  - DESCARTES hepatoblasts → *liver development*: BAAT and HMGCS2 lost the term.
  - KEGG Alzheimer's → *neuron apoptotic process*: GAPDH lost the term.
  - WP non-alcoholic fatty liver disease → *lipid metabolic process*: TNFRSF1A lost the term.
- `aigr_full` adds back 1: TRAVAGLINI airway smooth muscle → *smooth muscle contraction*, through new ACTA2/ACTG2 annotations.
- The 218 recovered under `aigr_prune` equals the lowest of the 20 null runs; only 2 of 20 nulls do as badly. That is a small loss, roughly equal to random removal of the same number of rows. The mean −log10 p of CORE terms falls by 0.12, about the same as the null.
- **The curated gold cannot show a precision gain.** Its `nonspecific` terms are hit only 3 times in any variant, so it cannot score noise reduction.

### 2. Pruning removes more hits than random removal, and the hits it removes are uninformative

The no-gold measure is: of the significant hits a variant loses relative to
`all`, how specific were they? Specificity is the mean information content
(IC) of the term, from baseline term sizes. "General" means the term annotates
more than 1,000 genes.

| arm | `aigr_prune` hits lost (net) | null hits lost (range over 20) | lost-hit mean IC: prune vs null | lost hits general: prune vs null |
|---|---|---|---|---|
| benchmark | 853 | 300–432 | **6.49** vs 6.76–7.68 | **32%** vs 17–29% |
| chem | 339 | 130–191 | **4.95** vs 5.02–5.67 | **61%** vs 46–60% |
| onco | 129 | 70–107 | 4.35 vs 4.26–4.96 | 72% vs 61–81% |

- **Benchmark and chemical exposures:** AIGR pruning removes about twice as many significant hits as random removal of the same rows. Its lost hits are more general than in every one of the 20 null runs. In chem, 102 sets end up with fewer terms and 28 with more.
- **Oncogenic (C6):** more hits are lost than in any null run, but their specificity falls inside the null range. These sets have the lowest reviewed-gene coverage (12%).
- **What is lost:** broad terms, e.g. *response to chemical*, *cellular response to chemical stimulus*, *response to endogenous stimulus*, *regulation of gene expression*, *regulation of catabolic process*, *protein-containing complex binding* and *membrane-bounded organelle*.
- **Sets whose result changes completely:**
  - One set, MASRI_…TAMOXIFEN_UP, loses its only hits (two broad developmental terms).
  - Two fenretinide sets (FERRARI_… and APPIERTO_RESPONSE_TO_FENRETINIDE_UP) gain hits they did not have. Pruning shrinks term sizes, which can sharpen a real signal.

### 3. Dropping non-core annotations hurts

Under `aigr_core`, recall_core falls from 0.565 to 0.535. Per set, 14 sets get worse and 2 get better (sign test p = 0.004). The mean −log10 p of CORE terms falls by 1.27. KEEP_AS_NON_CORE annotations carry real biology for enrichment, so a "core functions only" GOA would be a worse enrichment resource.

### 4. Scaled up: 16,622 MSigDB sets, split by review coverage

Every human MSigDB collection on MyGeneset except C5 (C5 is built from GO).
Sets have 10–1,000 genes. No gold: these are scored on hit counts and on how
specific the removed hits are, against the 20 random-removal nulls.

Breadth collections use a **per-set** threshold: raw p < 0.05 / 21,062 GO
terms (2.4e-6), the same for every collection and variant. genesets-rs's own
Bonferroni is matrix-wide (sets × terms), which would judge a 5,000-set
collection about 40× more strictly than the 128-set benchmark. The curated
benchmark in section 1 keeps the matrix-wide correction, so it stays
comparable with `evals/iba_vs_benchmark`. Outputs: `breadth_summary.tsv` and
`breadth_by_exposure.tsv`. The run takes about 80 minutes on 4 cores.

**The effect grows with review coverage.** Sets are pooled across all
collections and split by the share of their genes that AIGR has reviewed:

| reviewed share of set | sets | hits removed by `aigr_prune` (net) | by random removal (20 runs) | lost-hit mean IC: prune vs null | lost hits general: prune vs null |
|---|---|---|---|---|---|
| < 10% | 3,398 | 852 (0.25%) | 454–827 (0.13–0.24%) | **5.56** vs 5.85–6.26 | **51%** vs 38–46% |
| 10–25% | 11,374 | 37,148 (2.2%) | 14,605–18,313 (0.9–1.1%) | **5.67** vs 5.82–6.17 | **47%** vs 40–47% |
| ≥ 25% | 1,810 | 14,257 (**4.9%**) | 5,107–6,841 (1.7–2.3%) | **6.59** vs 7.07–7.63 | **37%** vs 22–30% |

Where AIGR has reviewed at least a quarter of a set's genes, pruning removes
about 5% of significant hits, 2–3× what random removal of the same rows does.
The removed hits are more general than in every null run. Below 10% coverage
the effect almost disappears. This supports the conclusion that the overall
effect is limited by coverage, not by the quality of the edits.

**By collection** (`aigr_prune` net hits lost vs the null range; specificity of the lost hits vs the null):

| collection | sets | median reviewed share | hits lost: prune vs null | lost-hit mean IC: prune vs null | lost hits general: prune vs null |
|---|---|---|---|---|---|
| Reactome | 1,497 | 16% | 5,658 vs 2,529–3,321 | **6.79** vs 7.38–7.88 | **34%** vs 18–26% |
| C2 pathways (BioCarta, PID, WikiPathways, KEGG) | 1,863 | 22% | 12,668 vs 4,497–6,413 | **6.81** vs 7.17–7.96 | **32%** vs 15–27% |
| C8 cell type | 829 | 11% | 1,955 vs 774–1,217 | **6.24** vs 6.56–6.94 | **35%** vs 25–32% |
| C3 TF/miRNA targets | 3,454 | 12% | 5,341 vs 2,183–3,133 | **5.12** vs 5.31–5.63 | **54%** vs 42–51% |
| C2 CGP (non-chemical) | 2,148 | 15% | 5,784 vs 2,223–3,045 | **5.99** vs 6.05–6.46 | 40% vs 31–40% |
| Hallmark | 50 | 20% | 444 vs 143–245 | 6.88 vs 6.73–7.70 | 27% vs 16–28% |
| benchmark (per-set threshold) | 128 | 17% | 1,269 vs 456–629 | 7.12 vs 7.05–7.80 | 24% vs 15–25% |
| C2 chemical exposure | 270 | 16% | 656 vs 216–361 | 5.69 vs 5.56–6.37 | 49% vs 36–50% |
| C4 cancer modules | 1,006 | 14% | 2,959 vs 941–1,327 | 6.16 vs 6.07–6.57 | 39% vs 32–42% |
| C7 immunologic | 5,148 | 14% | 15,322 vs 5,091–6,511 | 5.07 vs 5.04–5.37 | 59% vs 55–64% |
| C6 oncogenic | 189 | 12% | 201 vs 102–209 | 5.26 vs 5.04–5.72 | 52% vs 43–59% |

- **Volume:** in 10 of 11 collections pruning removes 2–3× more hits than random removal. C6, the smallest and least reviewed, is the exception.
- **Specificity:** the removed hits are more general than in every null run in Reactome, the C2 pathway databases, C8 and C3. C2 CGP passes on IC only. Hallmark, the benchmark and C2 chemical exposure sit at the edge of the null range on the per-set threshold (the chemical-exposure result in section 2 used the matrix-wide threshold). C4, C7 and C6 fall inside the null range.
- **Reactome caveat:** GOA carries Reactome-derived annotations, so the Reactome collection is not fully independent of GOA.
- **`aigr_full`** loses slightly more hits than `aigr_prune`, and the hits it loses are more specific. MODIFY swaps general terms for specific ones, which moves significance down the ontology.
- **`aigr_core`** removes 15–20% of all hits. Combined with its loss of CORE recall (section 3), that is too aggressive.

### With protein binding kept (`--drop-direct-terms ""`)

The first run of this eval kept `protein binding`. Dropping it matters a great deal for GOA and very little for the conclusions:

- It removes 278k of 905k NOT-filtered GAF rows (31%), but only 14k (gene, term) pairs. Most rows are IPI, many rows per pair.
- **About 85% of AIGR's removal rows were `protein binding`:** 10,311 of 11,421 REMOVE rows and 20,434 of 24,741 MARK_AS_OVER_ANNOTATED rows. Few of those deleted a pair, because a descendant binding term or an unreviewed row kept `protein binding` alive.
- 12 chem sets, 3 onco sets and 1 benchmark set (HP_SPINAL_MUSCULAR_ATROPHY) had `protein binding` (or `binding`) as their only significant result. Those sets now have none.
- **The AIGR effect does not depend on it.** With `protein binding` kept, pruning lost 912 benchmark hits against 368–507 for random removal, and 374 chem hits against 167–272. The specificity gap to the null was the same.
- CORE recall with `protein binding` kept was 221/394 for `all`, 219 for `aigr_prune` and 210 for `aigr_core`, over 126 scored sets.

## Interpretation

Applying AIGR's over-annotation removals to GOA makes enrichment results
**cleaner at little cost to recall of curated biology**. AIGR selectively strips
uninformative hits (broad response, regulation and compartment terms), well
beyond what random removal of the same number of annotations does. It does
**not** recover more CORE biology, and it loses about as much as random
removal does (3 of 391 terms). The effect is modest because only about 15% of
the genes in a typical external set have been reviewed, and because most
removed rows are backed by other rows asserting the same term.

## Review queue (not auto-applied)

`curated_core_term_changes.tsv` lists CORE terms lost through AIGR edits.
Check each against the review before deciding whether the review or the
enrichment expectation is wrong:

- BAAT and HMGCS2 → *liver development*. This is DESCARTES hepatoblasts; it is lost through MARK_AS_OVER_ANNOTATED.
- GAPDH → *neuron apoptotic process* (KEGG Alzheimer's).
- TNFRSF1A → *lipid metabolic process* (WP non-alcoholic fatty liver disease).
- *Lipid kinase activity* (`GO:0001727`) drops out of 12 KEGG cancer and diabetes sets. It is not a CORE term, and the PI3K biology is not lost. Take KEGG_GLIOMA:
  - In baseline the term rests on 6 members: PIK3CA, PIK3CB, PIK3CD, PIK3CG, PIK3R1 and PIK3R3 (p_adj 2.9e-3).
  - AIGR REMOVEs PIK3R1 (p85α) → *phosphatidylinositol kinase activity* (ISS). p85α is the regulatory subunit, not the kinase. With 5 members left the term misses the matrix-wide Bonferroni cutoff.
  - Under `aigr_prune` the specific signal is unchanged: *PI3K complex, class IA* (p_adj 1e-15), *PI3K/AKT signal transduction* (1e-19), *1-phosphatidylinositol-3-kinase activity* and *…regulator activity* all stay significant.
  - Follow-up: PIK3R3, also a regulatory subunit, still carries *1-phosphatidylinositol-3-kinase activity* (TAS, PMID:9524259, 2003). It has no AIGR review yet.

## Guardrail

As in `evals/iba_vs_benchmark`, the gold is never edited to match a variant.
A CORE term that an AIGR edit removes is either a review error or an
expectation error. Each one is a curator review item.
