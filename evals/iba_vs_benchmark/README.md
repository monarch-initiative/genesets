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

Per variant, over the CORE terms (`category` in `core_process`/`core_component`):

- **`recall_supported`** — recovered / total of `annotation_supported` core terms.
  The *fair* target: the curator asserts the genes carry these, so enrichment
  should recover them.
- **`gap_recovered`** — count of `annotation_gap` core terms recovered. The
  curator marked these too shallow in current GOA; recovery here means either a
  method surfaced it anyway or the gold label is conservative (a re-label
  candidate annotation_gap -> annotation_supported).
- **`unique_vs_baseline`** — supported-core terms a variant recovers that the
  `all` baseline does not.

## Headline result (2026, 98 sets, GOA `goa_human` current)

| variant | recall_supported | gap_recovered | unique_vs_all |
|---|---|---|---|
| all | 0.707 (145/205) | 0/40 | – |
| no_contributes_to | 0.702 (144/205) | 0/40 | 0 |
| iba_iea | 0.566 (116/205) | 0/40 | 0 |
| iba | 0.473 (97/205) | 0/40 | 2 |

(Post-fix numbers — see finding 3.)

1. **IBA carries ~2/3 of the core biology full GOA does** (0.47 vs 0.69); IEA
   recovers much of the difference (iba_iea 0.56).
2. **IBA is nearly a strict subset of all-GOA — it does not fill experimental
   gaps here.** Restricting to IBA loses 46 supported-core terms and uniquely
   recovers only 2, both conserved-housekeeping cellular components (ribosome
   `GO:0005840`, nucleolus `GO:0005730`).
3. **The eval audited and corrected the gold.** The first run found 9
   `annotation_gap` core terms that standard all-GOA enrichment *did* recover
   (e.g. SCHUHMACHER_MYC -> rRNA processing, TRAVAGLINI_CILIATED -> axoneme
   assembly, KEGG_RCC -> positive regulation of angiogenesis). Those were
   over-pessimistic gap predictions and were re-labeled `annotation_supported`;
   the table above is post-fix, with `gap_recovered` now 0/40 across every
   variant — the gold's `annotation_gap` set is empirically self-consistent (no
   remaining gap term is recoverable by any evidence variant).

Calibration: even all-GOA recovers only ~71% of `annotation_supported` core under
Bonferroni, so that label means "the genes carry it", not "it always reaches
genome-wide significance".

Generated tables (`/tmp/iba_eval/`) are not committed; rerun the pipeline to
regenerate. The scorer is `scripts/score_method_vs_benchmark.py`.
