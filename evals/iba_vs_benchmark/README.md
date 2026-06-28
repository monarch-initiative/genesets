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
   `_id`/`source=msigdb` hit. (The 31 `LIT:` sets store markers in prose, not
   membership, so they are out of scope here; 11 MSigDB sets — mostly newer C8
   single-cell clusters — are absent from mygeneset, leaving **98** evaluable.)

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

## Headline result (2026, 98 sets, GOA `goa_human` current)

| variant | recall_core | recall_supported | gap_recovered (disagreements) |
|---|---|---|---|
| all | 0.600 (180/300) | 0.694 | 9 |
| no_contributes_to | 0.597 (179/300) | 0.689 | 9 |
| iba_iea | 0.480 (144/300) | 0.561 | 6 |
| iba | 0.383 (115/300) | 0.469 | 5 |

1. **IBA carries ~2/3 of the core biology full GOA does** (recall_core 0.38 vs
   0.60); IEA recovers much of the difference (iba_iea 0.48).
2. **IBA is nearly a strict subset of all-GOA — it does not fill experimental
   gaps here.** Restricting to IBA loses 46 supported-core terms and uniquely
   recovers only 2, both conserved-housekeeping cellular components (ribosome
   `GO:0005840`, nucleolus `GO:0005730`).
3. **The eval surfaces a review queue — it does not edit the gold.** 9
   `annotation_gap` core terms were recovered by standard all-GOA enrichment
   (e.g. SCHUHMACHER_MYC -> rRNA processing, TRAVAGLINI_CILIATED -> axoneme
   assembly, KEGG_RCC -> positive regulation of angiogenesis). These are
   *disagreements* between the curator's gap prediction and a tool run, queued
   for deliberate curator review against GOA facts — not auto-relabeled.

Calibration: even all-GOA recovers only ~69% of `annotation_supported` core
under Bonferroni, so that label means "the genes carry it", not "it always
reaches genome-wide significance".

## Confirmatory vs mechanistic (the insight split)

With the corpus-wide `insight` tags, recall splits by whether a term is
`confirmatory` (restates the set's construction) or `mechanistic` (a non-obvious
process — a genuine enrichment insight). Over the evaluable sets:

| variant | recall_confirm (n=505) | recall_mechan (n=44) |
|---|---|---|
| all | 0.626 | 0.295 |
| no_contributes_to | 0.624 | 0.295 |
| iba_iea | 0.483 | 0.205 |
| iba | 0.364 | 0.227 |

**Mechanistic insight is ~2x harder to recover than confirmatory biology** —
even all-GOA recovers only ~30% of mechanistic terms vs ~63% of confirmatory
ones. Standard enrichment surfaces the obvious and largely misses the
non-obvious convergent mechanisms the curators flagged (often the
`annotation_gap` ones). The `mechanistic` denominator is small (44 evaluable of
the corpus's 50, since several live in `LIT:` sets with no fetched membership),
so those numbers are noisier — growing mechanistic-rich, agnostically-derived
sets is the way to sharpen this measure.

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
