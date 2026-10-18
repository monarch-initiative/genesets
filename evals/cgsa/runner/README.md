# cGSA local-runner adapter

cGSA has **no API**. It runs only as a local Python app that needs an Azure
OpenAI key. This adapter does the two halves we can own deterministically and
documents the manual run in between.

```
prepare   our GMT + sets TSV      ->  Data/<name>.xlsx  (cGSA input schema)
(manual)  run cGSA with your key   ->  cGSA_Results/<id>/*.json
collect   cGSA_Results/<id>/       ->  predictions TSV (result_id, pathway, score)
```

## 1. Prepare the input

```bash
python evals/cgsa/runner/run_cgsa.py prepare \
  --gmt  evals/cgsa/data/case_study_members.gmt \
  --sets evals/cgsa/data/case_study_sets.tsv \
  --out  /tmp/cgsa_input.xlsx
```

Writes an `.xlsx` with the sheet **Sheet1** and cGSA's expected columns:
`PMC, Topic, Database, Context, Gene List, Species` (Database is one of cGSA's
Enrichr-library keys — GOBP/GOMF/GOCC/KEGG/MSigDB/Reactome/Enrichment/Pathway/
Disease/Family; Species is 9606 or 10090). Works with either the case-study
files or `benchmark_members.gmt` + `benchmark_sets.tsv`.

## 2. Run cGSA (manual — needs the repo + an Azure OpenAI key)

```bash
git clone https://github.com/ncbi-nlp/cGSA && cd cGSA
# set openai.api_key / api_base / api_version in:
#   evaluation.py, exploration.py, confidence.py
cp /tmp/cgsa_input.xlsx Data/
# point PathDis.py's read_excel(...) at Data/cgsa_input.xlsx (it defaults to demo.xlsx)
python PathDis.py
# -> cGSA_Results/<id>/{Communities,Tradition_GSEA,Context_GSEA,
#                        Summarized_Results,Summarized_Results_Refine_with_Scores}.json
```

cGSA also needs `p_threshold` (0.001) and `ptw_threshold` (6.0) as set in
PathDis.py, plus its Python 3.11 dependency set (see the upstream README).

## 3. Collect the output

```bash
python evals/cgsa/runner/run_cgsa.py collect \
  --results-dir <cGSA>/cGSA_Results \
  --out /tmp/cgsa_preds.tsv
```

Reads each `Summarized_Results_Refine_with_Scores.json` and emits
`result_id, cgsa_pathway, relevance_score` — the same schema as
`evals/cgsa/data/cgsa_predictions.tsv` (minus the expert category, which only the
authors' blinded annotation supplies). The parser is defensive about cGSA's
evolving JSON shape.

## 4. Score

Ground the collected free-text pathways to GO IDs with
`scripts/ground_cgsa_to_go.py`, then feed the resulting `query_id/target_id`
file to `scripts/score_method_vs_benchmark.py` (or compare against the grounded
benchmark ground truth). See `../README.md`.
