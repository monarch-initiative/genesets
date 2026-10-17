# GSEA methods against the curated benchmark

Scores `genesets-rs gsea` configurations against the curator CORE terms in
`curation/genesets/*.yaml`, the same gold that `evals/iba_vs_benchmark` uses
for Fisher. The gold's sets carry no scores, so each set gets synthetic ones:
members draw from N(shift, 1), at 1.5 and 3.0 SD, and every other annotated
gene from N(0, 1). Fisher on the set itself (the existing benchmark pipeline)
is the reference; Fisher on the top 5% of the synthetic ranking is the
thresholded analogue.

This measures how the statistics and the graph-aware algorithms trade recall of
curator-level biology against the number of terms they report. It does not
measure performance on real expression data. The gold is only read, never
refit (see the guardrail in `evals/iba_vs_benchmark/README.md`).

```bash
python3 scripts/prepare_go_eval.py --out-dir /tmp/go_eval --variants all
cargo build --release
uv run --project python/genesets-workflows --extra curation \
  python evals/gsea_benchmark/run_benchmark_eval.py --go-dir /tmp/go_eval \
  --work-dir /tmp/gsea_benchmark --out evals/gsea_benchmark/RESULTS.md
```

Membership is fetched from mygeneset.info (MSigDB sets) and
`curation/genesets/lit_members.gmt` (literature sets) into
`<work-dir>/queries.gmt`, then reused. Results are in [RESULTS.md](RESULTS.md).

## Findings (October 2026 run; GO 2026-04-27, GOA 2026-05-21)

From [RESULTS.md](RESULTS.md), 126 scored sets with 394 CORE terms:

- **Ranked methods recover far less than Fisher on the exact set** (recall
  0.561). At a 3 SD shift KS classic is the best ranked method (0.226),
  followed by Fisher on the top 5% (0.152) and Wilcoxon classic (0.102); at
  1.5 SD all are below 0.08. Exact membership is an oracle no ranking matches.
- **KS beats Wilcoxon here**, the reverse of the planted eval's top-rank
  result: a GO term holds a set's genes plus many others, and KS's weighted
  running sum picks up a subset at the top that Wilcoxon's whole-term rank
  shift dilutes.
- **elim and weight01 roughly halve recall of curator CORE terms** (KS at
  3 SD: 0.226 to 0.114 and 0.147) while cutting the terms reported per set
  (mean 60 to 19 and 37). Curators often place CORE biology at levels these
  algorithms strip as redundant. weight01 keeps more than elim throughout.
- **Curated real rankings.** Entries with a `ranking` block are also scored on
  their real scores, in a separate table. So far that is one entry,
  `LIT:SCHIZOPHRENIA_SCHEMA_EXOME` (rare-variant burden, unsigned). Its CORE
  terms were seeded by Wilcoxon classic, so that row's recall of 1.0 is not an
  independent test. The informative rows: KS recovers none of them (the signal
  is diffuse, not a concentrated top), Wilcoxon weight01 keeps half, and
  Wilcoxon elim none. One entry is an illustration, not a result; the table
  becomes a benchmark as ranked entries accumulate.
