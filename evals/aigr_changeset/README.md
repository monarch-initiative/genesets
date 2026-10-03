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

**Change-set coverage.** 105,086 reviewed (non-NOT) annotation rows; 98,195
(93%) match a row in current GOA. 1,470 NEW rows. Reviewed genes make up a
median **17%** of the members of a benchmark set (16% chem, 12% onco).

Row removals mostly do not remove pairs. 11,421 REMOVE rows delete only 1,122
(gene, term) pairs, and REMOVE + MARK_AS_OVER_ANNOTATED (36,162 rows) delete
5,187 pairs. In most cases another, unreviewed or accepted, row still asserts
the same pair.

### 1. Recall of the curated biology is unchanged by pruning

| variant | recall_core (126 sets) | recall_confirm | recall_mechan | significant terms |
|---|---|---|---|---|
| all | 0.561 (221/394) | 0.569 | 0.333 | 29,032 |
| aigr_remove | 0.558 (220) | 0.567 | 0.333 | 28,937 |
| **aigr_prune** | **0.556 (219)** | 0.566 | 0.333 | **28,120** |
| aigr_full | 0.558 (220) | 0.564 | 0.321 | 28,100 |
| aigr_core | 0.533 (210) | 0.531 | 0.296 | 23,670 |
| null_prune (mean of 20) | 0.558 (219.8; range 218–221) | 0.566 | 0.333 | 28,577 (28,525–28,664) |

- `aigr_prune` loses 2 CORE terms:
  - DESCARTES hepatoblasts → *liver development*: BAAT and HMGCS2 lost the term.
  - KEGG Alzheimer's → *neuron apoptotic process*: GAPDH lost the term.
- `aigr_full` adds back 1: TRAVAGLINI airway smooth muscle → *smooth muscle contraction*, through new ACTA2/ACTG2 annotations.
- This is indistinguishable from random removal. The mean −log10 p of CORE terms falls slightly (−0.12), similar to the null. In sets with ≥25% reviewed members it falls a little more than in 19/20 nulls.
- **The curated gold has no way to show a precision gain.** Its `nonspecific` terms are hit only 3 times in any variant, so it cannot score noise reduction.

### 2. Pruning removes more hits than random removal, and the hits it removes are uninformative

The no-gold measure is: of the significant hits a variant loses relative to
`all`, how specific were they? Specificity is the mean information content
(IC) of the term, from baseline term sizes. "General" means the term annotates
more than 1,000 genes.

| arm | `aigr_prune` hits lost | null hits lost (range over 20) | lost-hit mean IC: prune vs null | lost hits general: prune vs null |
|---|---|---|---|---|
| benchmark | 912 net (28,120 sig) | 368–507 net | **6.40** vs 6.72–7.36 | **34%** vs 19–27% |
| chem | 374 net (12,623 sig) | 167–272 net | **4.77** vs 4.82–5.71 | **63%** vs 46–61% |
| onco | 120 net (11,794 sig) | 56–138 net | 4.39 vs 4.38–4.98 | 72% vs 63–79% |

- **Benchmark:** AIGR pruning removes about twice as many significant hits as random removal of the same rows. Its lost hits are more general than in every one of the 20 null runs.
- **Chemical exposures:** the same pattern. Removed hits are more general than in all 20 nulls, and 113 sets have fewer terms against 32 with more.
- **Oncogenic (C6):** the effect falls inside the null range. These sets have the lowest reviewed-gene coverage (12%).
- **What is lost:** the most-lost hits are `protein binding` (14 chem sets), `binding`, `extracellular exosome`, *regulation of gene expression*, *response to chemical* and similar terms.
- **Sets left with no hits:** six sets lose their only hits. Five had only `protein binding` (or `binding`) as their entire result: HP_SPINAL_MUSCULAR_ATROPHY, CHEOK_…MERCAPTOPURINE_DN, STAMBOLSKY_…VITAMIN_D3_DN, KANG_FLUOROURACIL_RESISTANCE_UP and CAFFAREL_RESPONSE_TO_THC_UP. The sixth, MASRI_…TAMOXIFEN_UP, loses two broad developmental terms.

### 3. Dropping non-core annotations hurts

Under `aigr_core`, recall_core falls from 0.561 to 0.533. Per set, 13 sets get worse and 2 get better (sign test p = 0.007). The mean −log10 p of CORE terms falls by 1.26. KEEP_AS_NON_CORE annotations carry real biology for enrichment, so a "core functions only" GOA would be a worse enrichment resource.

## Interpretation

Applying AIGR's over-annotation removals to GOA makes enrichment results
**cleaner without costing recall of curated biology**. AIGR selectively strips
the uninformative hits (`protein binding`, `exosome`, broad regulation terms),
and does so more than random removal of the same number of annotations. It
does **not** recover more CORE biology. The effect is modest because only
about 15% of the genes in a typical external set have been reviewed, and
because most removed rows are backed by other rows asserting the same term.

## Review queue (not auto-applied)

`curated_core_term_changes.tsv` lists CORE terms lost through AIGR edits.
Check each against the review before deciding whether the review or the
enrichment expectation is wrong:

- BAAT and HMGCS2 → *liver development*. This is DESCARTES hepatoblasts; it is lost through MARK_AS_OVER_ANNOTATED.
- GAPDH → *neuron apoptotic process* (KEGG Alzheimer's).
- *Lipid kinase activity* (`GO:0001727`) drops out of 12 KEGG cancer and diabetes sets. It is not a CORE term, and the PI3K biology is not lost. Take KEGG_GLIOMA:
  - In baseline the term rests on 6 members: PIK3CA, PIK3CB, PIK3CD, PIK3CG, PIK3R1 and PIK3R3 (p_adj 2.9e-3).
  - AIGR REMOVEs PIK3R1 (p85α) → *phosphatidylinositol kinase activity* (ISS). p85α is the regulatory subunit, not the kinase. With 5 members left the term misses the matrix-wide Bonferroni cutoff.
  - Under `aigr_prune` the specific signal is unchanged: *PI3K complex, class IA* (p_adj 1e-15), *PI3K/AKT signal transduction* (1e-19), *1-phosphatidylinositol-3-kinase activity* and *…regulator activity* all stay significant.
  - Follow-up: PIK3R3, also a regulatory subunit, still carries *1-phosphatidylinositol-3-kinase activity* (TAS, PMID:9524259, 2003). It has no AIGR review yet.

## Guardrail

As in `evals/iba_vs_benchmark`, the gold is never edited to match a variant.
A CORE term that an AIGR edit removes is either a review error or an
expectation error. Each one is a curator review item.
