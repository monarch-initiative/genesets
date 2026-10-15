# How should enrichment treat "regulation of X"?

**Question:** GO keeps "regulation of X" terms separate from X. With the
standard is_a + part_of closure, a gene annotated to *positive regulation of
angiogenesis* reaches *regulation of angiogenesis* but never *angiogenesis*.
Does enrichment recover more of the curated biology if we (a) drop the
regulation terms, (b) follow `regulates` edges in the closure, or (c) collapse
regulation annotations onto the regulated process?

Companion to the annotation-side audit in ai-gene-review
(`projects/REGULATION_CONSISTENCY.md`), which asks how much of GOA's
regulation annotation is internally inconsistent.

## Variants

All four share NOT-filtered `goa_human`, one all-GOA background, Bonferroni
`max_p_adjust 0.05`, `min_overlap 2`. Only the closure and/or the gene-term
table differ.

| variant | closure | regulation annotations | regulation terms tested |
|---|---|---|---|
| `current` | is_a + part_of | kept | yes |
| `drop_regulation` | is_a + part_of | removed | no |
| `regulates_closure` | is_a + part_of + regulates (+/−/unsigned) | kept | yes |
| `collapse_to_target` | is_a + part_of | re-pointed to the regulated X | no |

`current` is what `prepare_go_eval.py` does and what the other benchmark evals
use. Sign is already pooled in `current`: *positive/negative regulation of X*
are `is_a` *regulation of X*, so a sign-agnostic regulation term needs no new
variant.

**Defining a regulation term.** A BP term with a direct `regulates`,
`positively_regulates` or `negatively_regulates` edge, plus regulation-named
`is_a` descendants ("regulation of X by Y"). The name restriction matters: GO
classifies *signal transduction* — and so every signaling pathway — `is_a`
*regulation of cellular process*. Following `is_a` without it labels all
signaling as regulation and makes `drop_regulation` look far worse than it is
(recall_core 0.439 rather than 0.513 in a first run). MF regulator edges
(*enzyme inhibitor activity* `negatively_regulates` *catalytic activity*) are
out of scope.

## Run

```bash
E=/tmp/reg_eval
python3 scripts/fetch_benchmark_queries.py --out $E/queries.gmt
python3 scripts/prepare_regulation_variants.py --out-dir $E   # downloads full go.obo + goa_human
for v in current drop_regulation regulates_closure collapse_to_target; do
  genesets-rs run $E/$v/config.yaml
done
uv run --project python/genesets-workflows --extra curation \
  python scripts/score_regulation_variants.py --eval-dir $E --obo $E/downloads/go.obo \
  --out $E/scores.tsv --changes-out $E/changes.tsv
```

`prepare_regulation_variants.py` needs the full `go.obo`; go-basic drops the
regulates relations.

## Results (GO 2026-07-26, goa_human 2026-05-21; 126 scored sets, 394 CORE terms)

[`scores.tsv`](scores.tsv), per-set gains/losses in [`changes.tsv`](changes.tsv).
`current` reproduces the `iba_vs_benchmark` `all` baseline exactly
(recall_core 0.561, 221/394).

| variant | hits | reg share of hits | recall_core | CORE participant terms (366) | CORE regulation terms (28) | lenient | curated-negative hits |
|---|---|---|---|---|---|---|---|
| current | 29,068 | 31% | 0.561 | 0.552 | 0.679 | 0.569 | 3 |
| drop_regulation | 20,224 | 0% | 0.513 | 0.552 | 0 | 0.538 | 3 |
| regulates_closure | 33,046 | 27% | **0.617** | **0.612** | 0.679 | **0.627** | 4 |
| collapse_to_target | 24,292 | 0% | 0.566 | 0.609 | 0 | 0.624 | 4 |

*lenient*: a gold *regulation of X* term also counts as recovered when the
variant hits X. That is the fair score for variants that cannot report
regulation terms.

1. **Dropping regulation costs nothing on the participant biology.** Participant
   recall is unchanged (0.552), and 30% fewer terms reach significance. All 19
   lost CORE terms are regulation terms: *regulation of cell population
   proliferation* in 12 KEGG cancer sets, *regulation of transcription by RNA
   polymerase II* (4), and one each of *negative regulation of cell population
   proliferation*, *positive regulation of angiogenesis* and *negative
   regulation of insulin receptor signaling pathway*.
2. **Regulates closure is the only variant that improves strict recall**
   (+0.056, +22 CORE terms net). The gains are participant terms that
   regulators push over threshold: *T cell mediated cytotoxicity* in five
   immune KEGG sets, *insulin secretion* and *cellular response to glucose* in
   beta cells, *EMT* in the hallmark EMT set, *ferroptosis* in the
   ferroptosis-suppressor panel, *amyloid precursor protein catabolism* in the
   Alzheimer GWAS genes. It costs 14% more output and one more curated-negative
   hit.
3. **Collapsing gets most of the participant gain with a smaller output.**
   `collapse_to_target` matches the regulates closure on participant recall
   (0.609 vs 0.612) and on lenient recall (0.624 vs 0.627), with 27% fewer hits
   than the closure and 16% fewer than `current`.
4. **Which variant is best depends on what you want reported.** If a regulation
   term is an acceptable answer, use the regulates closure. If you want one
   answer per process and a shorter list, collapse. Dropping is never better
   than collapsing on these metrics.

## Recuration: pooled, blinded judgments

The gold above leans towards `current`. 880 of its 1,063 associations were
seeded from enrichment results, presumably under the default is_a/part_of
closure, so terms reachable only through `regulates` were less likely to be
proposed. To remove that bias without refitting the gold to any one variant,
the contested region was recurated the way IR benchmarks fix pooling bias.

1. **Pool symmetrically.** For each set, take the top 3 hits by p-value from
   *each* variant, restricted to terms the variants disagree on (significant
   in some but not all) with at most 1,000 annotated genes, minus terms
   already curated. This gives 884 candidates over 115 sets
   ([`recuration/build_pool.py`](recuration/build_pool.py),
   [`recuration/pool.json`](recuration/pool.json)).
2. **Judge blind.** Each candidate was shown with the set's full membership,
   the member genes annotated to it via is_a/part_of, and the members reaching
   it only via `regulates`, but **not** which variant proposed it. Judges
   assigned `category`, `recovery_status` and `insight` under the
   `curate-geneset` rules, negatives included. The 884 judgments are in
   [`recuration/judgments.yaml`](recuration/judgments.yaml): 83 core_process,
   4 core_component, 297 supporting, 185 marker-driven, 232 nonspecific, 83
   false_association.
3. **Rescore** with the judgments added to a copy of the gold
   (`scripts/merge_pooled_judgments.py --genesets-dir <copy>`). The judgments
   are **not yet merged into `curation/genesets/`**; that is a curator decision.

### Recall against the recurated gold (481 CORE terms, 60 of them regulation terms)

[`scores_recurated.tsv`](scores_recurated.tsv)

| variant | recall_core | participant CORE | regulation CORE | lenient | recall_confirm | curated-negative hits |
|---|---|---|---|---|---|---|
| current | 0.549 | 0.506 | 0.850 | 0.555 | 0.578 | 180 |
| drop_regulation | 0.459 | 0.525 | 0 | 0.516 | 0.408 | 118 |
| regulates_closure | **0.653** | 0.625 | 0.850 | **0.661** | **0.687** | 215 |
| collapse_to_target | 0.555 | **0.634** | 0 | **0.661** | 0.518 | 149 |

### Precision in the contested region

Of the judged candidates each variant hits ([`pooled_precision.tsv`](pooled_precision.tsv),
[`recuration/score_pooled_precision.py`](recuration/score_pooled_precision.py)):

| variant | judged hits | core | core + supporting | nonspecific / false | positive rate |
|---|---|---|---|---|---|
| current | 498 | 43 | 228 | 177 | 0.458 |
| drop_regulation | 264 | 19 | 94 | 115 | 0.356 |
| regulates_closure | 643 | 71 | 297 | 211 | **0.462** |
| collapse_to_target | 408 | 44 | 163 | 145 | 0.400 |

### What changed

1. **Regulation terms are not the noise.** In the contested region, judges
   called 50% of regulation-term candidates core or supporting (149/300),
   against 40% of participant-term candidates (235/584). The first-pass result
   that dropping regulation was free came from a gold that had rarely been
   offered regulation terms. Against the recurated gold, dropping costs 9
   points of recall and has the lowest contested-region precision.
2. **Regulates closure wins outright.** It has the best strict and lenient
   recall, and the most positive hits (297 vs 228 for `current`) at the same
   positive rate. Its extra output is as good as the output it already had,
   though in absolute terms it also adds the most curated negatives.
3. **Collapse is the compact alternative.** It ties regulates closure on
   participant and lenient recall with 27% less output, but it cannot report
   the 60 core regulation terms, and its contested-region precision is lower.

The recommendation therefore changes from "collapse or drop" to "regulates
closure, or collapse if you need one answer per process".

## Caveats

- **Pool depth.** Only the top 3 contested hits per variant per set were
  judged, so recall on pool-added terms reflects that depth. The pool is
  symmetric, so no variant gets more nominations than another.
- **The judges were LLM curators** working under the `curate-geneset` rules,
  one judge per candidate. There is no inter-annotator agreement yet. A
  spot-check by a human curator should precede merging into the gold.
- **Validation of the merged copy was partial.** LinkML schema validation
  passed, and all 884 GO ids were checked against `go.obo` 2026-07-26 for
  exact label and non-obsolescence. The full `curate validate` could not run
  in the authoring environment because its OAK step failed to download
  `mondo.db`. Run `just curate-validate` before merging.
- **The original gold is unchanged.** The first-pass table above still scores
  against it.
- **Curated-negative hits are a floor, not a false-positive rate.** Most hits
  outside the pool were never curated either way.
- **Regulates closure inherits annotation noise.** It also propagates the
  sign-conflicting and participant-misfiled regulation annotations that the
  ai-gene-review audit flags. A cleaner GOA would change the closure results
  more than the others.
- Regulation-of-regulation chains are followed to the first non-regulation
  process, as AmiGO's `regulates_closure` does. That transitivity is a
  convenience, not an RO entailment.
