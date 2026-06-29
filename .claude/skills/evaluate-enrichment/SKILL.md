---
name: evaluate-enrichment
description: >-
  Use when evaluating a GO enrichment method or a GOA evidence variant (all vs
  IBA vs IBA+IEA vs no-contributes_to) against the curated benchmark in this repo
  — the IBA/evidence-ablation eval. Covers building queries.gmt, prepare_go_eval,
  running genesets-rs, and scoring with the confirmatory/mechanistic split and the
  recovery_status diagnostics, plus the don't-refit-the-gold guardrail.
---

# Evaluating enrichment against the curated benchmark

This scores an enrichment run's GO terms against external ground truth: the
curator-validated CORE terms in `curation/genesets/*.yaml`. Unlike the
`evals/go_iba_impact_*` analyses (which compare annotation variants to *each
other*, no ground truth), this is a real precision/recall measurement.

**Read first:** `evals/iba_vs_benchmark/README.md` — the full pipeline, the
current headline numbers, and the guardrail. This skill is the runbook.

## Pipeline

All generated artifacts live under a scratch dir (e.g. `/tmp/iba_eval/`) and are
NOT committed; rerun to regenerate.

1. **Build `queries.gmt`** — the member genes of every evaluable set.
   - MSigDB sets: fetched from mygeneset.info (`_id` = set name) → a `queries.msigdb`
     base. See `scripts/fetch_mygeneset_query.py`.
   - `LIT:` sets: `curation/genesets/lit_members.gmt` (captured from primary sources).
   - `queries.gmt` = `queries.msigdb` + `lit_members.gmt`. Each line:
     `<gene_set_name>\t<desc>\tGENE1\tGENE2…`. The first field must equal the YAML
     `gene_set_name` (the scorer joins on it).

2. **Prepare annotation variants** — downloads GO + `goa_human.gaf.gz`, writes a
   per-variant evidence-filtered `gene_terms.tsv`, shared `terms.tsv`/`closure.tsv`,
   a common background, and a `genesets-rs` config per variant (Bonferroni, `max_p_adjust 0.05`):
   ```bash
   python3 scripts/prepare_go_eval.py --out-dir /tmp/iba_eval \
     --variants all,iba,iba_iea,no_contributes_to
   ```
   Place `queries.gmt` (step 1) at `/tmp/iba_eval/queries.gmt`.

3. **Run enrichment** per variant → `<variant>/results.tsv`:
   ```bash
   for v in all iba iba_iea no_contributes_to; do
     genesets-rs run /tmp/iba_eval/$v/config.yaml
   done
   ```

4. **Score vs the gold:**
   ```bash
   uv run --project python/genesets-workflows --extra curation \
     python scripts/score_method_vs_benchmark.py --eval-dir /tmp/iba_eval \
     --out /tmp/iba_eval/benchmark_scores.tsv
   ```

## Metrics (per variant)

- **`recall_core`** — recovered / total CORE terms (`category` in
  core_process/core_component). The headline; tool-independent of recovery_status.
- **`recall_confirm` (n=…)** / **`recall_mechan` (n=…)** — recall split by the
  `insight` tag. Mechanistic insight runs ~2× harder to recover than confirmatory;
  read `recall_mechan` alongside its (small) denominator.
- **`gap_recovered`** — count of `annotation_gap` CORE terms a variant recovered:
  a DISAGREEMENT between the curator's gap prediction and a tool run. A review
  item, NOT an automatic gold edit.

For per-set diagnostics (which terms a set recovered, by insight/recovery_status),
load the gold via `genesets_workflows.curation.model.load_interpretation` and the
hits from `<variant>/results.tsv` (`query_id` → `target_id`); see the inline
scripts in `evals/iba_vs_benchmark/README.md`.

## Adding new sets to the eval

When the **curate-geneset** skill adds sets, fold them in before re-scoring:
append MSigDB memberships to `queries.msigdb` and `LIT:` memberships are already
in `lit_members.gmt`; rebuild `queries.gmt = queries.msigdb + lit_members.gmt`;
re-run steps 3–4. The GO+GOA download (step 2) can be reused.

## The guardrail — the eval measures, it never refits the gold

`category` (biology) is the scored authority and stays fixed. `recovery_status`
and `insight` are curator judgments checked against GOA facts during curation —
**never** auto-updated to match whatever a method + GOA snapshot + p-threshold
recovers. Refitting them would make recall-vs-gold circular (the tool couldn't be
"wrong" because truth would be redefined to match it) and forfeit the gold's
independence across methods and GOA versions. `gap_recovered` disagreements are
inputs to deliberate curator review, not relabels. If you find yourself editing
`curation/genesets/*.yaml` to raise a recall number, stop — that is the failure
mode this guardrail exists to prevent.

Calibration note: even all-GOA recovers only ~56% of CORE terms under Bonferroni,
so `core` + `annotation_supported` means "the genes carry it", not "it always
reaches genome-wide significance".
