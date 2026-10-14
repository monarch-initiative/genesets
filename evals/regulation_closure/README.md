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

## Caveats

- **The gold leans towards `current`.** 880 of 1,063 curated associations were
  seeded from enrichment results (`seed_source: enrichment_recovered`). If that
  seeding used the standard is_a/part_of closure (the repo default), terms
  reachable only through `regulates` were less likely to be proposed, and the
  regulates-closure gains are an underestimate. Per the benchmark guardrail, the gold is not
  edited to suit any variant here.
- **Curated-negative hits are a floor, not a false-positive rate.** Most hits
  were never curated either way.
- **Regulates closure inherits annotation noise.** It also propagates the
  sign-conflicting and participant-misfiled regulation annotations that the
  ai-gene-review audit flags. A cleaner GOA would change the closure results
  more than the others.
- Regulation-of-regulation chains are followed to the first non-regulation
  process, as AmiGO's `regulates_closure` does. That transitivity is a
  convenience, not an RO entailment.
