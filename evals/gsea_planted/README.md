# Planted-signal GSEA evaluation

The curated benchmark interprets gene sets that carry no scores, so it cannot
on its own score ranked-list methods. This eval plants a known signal in a
known GO Biological Process term and checks how each `genesets-rs gsea`
configuration credits it, on the real GO graph and GOA human annotations.

- **specific**: half the genes of a 20-60 gene term shift up. A good method
  finds the term, ranks it first, and does not drag its ancestors along.
- **general**: half the genes of a 100-300 gene term with at least two
  children shift up. The signal belongs to the general term, so a method that
  pushes credit down to children loses it. This is the case elim is
  criticised for.
- **null**: no signal, to check calibration.

Methods: KS (GSEA) and Wilcoxon rank-sum, each with `classic`, `elim` and
`weight01`, plus Fisher ORA on the top 5% of the ranking as the thresholded
baseline. Design details are in the script's docstring; results are in
[RESULTS.md](RESULTS.md).

```bash
python3 scripts/prepare_go_eval.py --out-dir /tmp/go_eval --variants all
cargo build --release
python3 evals/gsea_planted/run_planted_eval.py --go-dir /tmp/go_eval \
  --work-dir /tmp/gsea_planted --out evals/gsea_planted/RESULTS.md
```

Generated rankings and outputs (about 1 GB) stay in `--work-dir` and are not
committed.

## Findings (October 2026 run; GO 2026-04-27, GOA 2026-05-21)

Numbers below are from [RESULTS.md](RESULTS.md): 100 specific, 100 general and
50 null replicates; 31,552 BP-annotated genes ranked, 4,425 BP terms tested.

- **The graph algorithms remove ancestor redundancy.** With a specific signal,
  classic KS flags 0.79 of the planted term's ancestors per run, classic
  Wilcoxon 0.39, and ORA 1.31. elim takes that to 0; weight01 to 0.18 (KS) and
  0.04 (Wilcoxon).
- **elim pays for it when the signal is general.** It finds the planted
  general term in 6% (KS) and 8% (Wilcoxon) of runs, against 66% and 85% for
  classic, with a median rank of 45-58: significant children take the genes
  and leave the parent empty. It also costs power on specific terms, whose own
  descendants cross the 0.01 cutoff (KS 0.78 to 0.44, Wilcoxon 0.51 to 0.18).
- **weight01 keeps most of a general signal with Wilcoxon.** It finds the
  general term in 56% of runs with a median rank of 1, while cutting
  significant ancestors from 0.85 to 0.02 and unrelated terms from 1.90 to
  0.26 per run. With KS it does much worse (16%): at the permutation floor
  parent and child p-values tie, and children usually have the larger |NES|.
- **KS versus Wilcoxon, classic.** KS finds more planted terms (specific 0.78
  vs 0.51) but rarely ranks them first (top1 0.05 vs 0.52): subsets that hold
  the shifted genes, descendants or overlapping terms, outrank the planted
  term on the weighted running sum.
- **Calibration.** With no signal, Wilcoxon and ORA report nothing; KS with
  the GSEA FDR reports 0.20 terms per null ranking at q <= 0.05, so the GSEA
  FDR is somewhat liberal here.
- **Thresholded ORA is competitive in this design.** A 2 SD shift on half the
  genes puts them in the top 5%, which suits a threshold. Do not read this as
  ORA beating GSEA on real data.
