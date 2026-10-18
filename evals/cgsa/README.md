# cGSA: benchmark import + evaluation adapters

Adapters and imported data for **cGSA** — *Knowledge-guided contextual gene set
analysis with large language models* (Wang, Day, Wei et al., **Bioinformatics**
2026, 42, btag214; code https://github.com/ncbi-nlp/cGSA). cGSA is an LLM
pipeline that clusters DEGs on a STRING PPI network (EdMot), runs Enrichr per
cluster, then uses GPT-4o to screen and summarize pathways against a stated
research **context**, emitting a short list of free-text pathway names with a
0–10 relevance score.

This directory brings cGSA's benchmark into the repo and lets us **evaluate
cGSA's approach** here, two ways.

## Is there an API? No.

cGSA has no hosted service or programmatic endpoint. It is a local Python app
(`python PathDis.py`) that requires an **Azure OpenAI key hard-coded into three
of its source files** plus live STRING / Enrichr / OLS / PubMed access. So
"evaluate their approach" means either (a) scoring the outputs the authors
already published, or (b) running their repo locally with your own key. Both are
supported below; only (a) runs unattended in this repo.

## Data (`data/`)

`data/build_cgsa_data.py` fetches the four supplementary workbooks from
`raw.githubusercontent.com/ncbi-nlp/cGSA` and flattens them to committed,
diff-friendly TSV/GMT. See `data/PROVENANCE.md`. Re-run to refresh.

| file | what |
|------|------|
| `benchmark_sets.tsv` | 102 benchmark DEG sets: contexts (objective/keywords/conditions), `database`, taxon, counts |
| `benchmark_ground_truth.tsv` | 1,671 paper-derived ground-truth pathway names (free text), per set |
| `benchmark_members.gmt` | the DEG membership (enrichment input) per set |
| `case_study_sets.tsv` / `case_study_members.gmt` | the 9 case-study DEG sets |
| `cgsa_ground_truth.tsv` / `cgsa_predictions.tsv` | the published results: 1,671 ground-truth + 1,482 cGSA pathways, each with the expert **functional-category** number |

Key structural facts: ground truth and cGSA outputs are both **free text** (not
GO IDs); the benchmark is 52% mouse; DEG sets average ~790 genes; `database`
tells you which Enrichr library each set was scored against (36 KEGG, 20 GOBP,
…). The two id spaces differ — the benchmark sheet keys sets by `PMC…`, the
results sheet by `PMID…` — so the published-outputs scorer is self-contained.

## Adapter 1 — score cGSA's published outputs (model-free)

`scripts/score_cgsa_published.py` reproduces the paper's **expert** HIT/ACC from
the functional-category annotation — a shared category between a ground-truth
pathway and a cGSA pathway is a match — needing no embedding model and no API:

```bash
python scripts/score_cgsa_published.py
```
```
result sets: 102   ground-truth pathways: 1671   cGSA pathways: 1482   off-context (cat 0): 126
variant  sets    ACC   HIT  harmonic_of_means  mean_per_set_harmonic
matched   102  0.837  0.86              0.849                  0.821
 strict   102  0.797  0.86              0.828                  0.798
paper Fig. 4a (expert): ACC 0.828  HIT 0.873  harmonic 0.850
```

`matched` reproduces the paper's headline **0.850**; `strict` keeps cGSA's
off-context (category-0) outputs in the precision denominator, the honest ACC
including irrelevant outputs. Counts (1671 / 1482) match the paper exactly.

## Adapter 2 — run cGSA locally on our gene sets (scaffold)

`runner/run_cgsa.py` owns the two deterministic halves and documents the manual
step in between (see `runner/README.md`):

```bash
# our GMT + sets TSV  ->  cGSA's Data/<name>.xlsx input schema
python evals/cgsa/runner/run_cgsa.py prepare \
  --gmt evals/cgsa/data/case_study_members.gmt \
  --sets evals/cgsa/data/case_study_sets.tsv --out /tmp/cgsa_input.xlsx
# (manual) clone ncbi-nlp/cGSA, set Azure keys, `python PathDis.py`
# cGSA_Results/<id>/  ->  predictions TSV
python evals/cgsa/runner/run_cgsa.py collect --results-dir <cGSA>/cGSA_Results --out preds.tsv
```

## Bridge — ground free-text pathways to GO IDs

`scripts/ground_cgsa_to_go.py` resolves free-text pathway names to GO IDs with
OLS exact-label lookup (cached to `data/go_grounding_cache.json` for
determinism). This turns cGSA output — or the benchmark's own ground truth —
into the `query_id<TAB>target_id` schema that `score_method_vs_benchmark.py`
already reads, so cGSA and `genesets-rs` become comparable on GO-ID terms. Only
GO-flavored names ground; KEGG/Reactome/MSigDB-signature names are reported
unmapped, never forced. Coverage on the 20 GOBP-database benchmark sets is
recorded in `data/benchmark_go_grounding.tsv`.

## Curated subset

Three DEG modules from one benchmark melanoma study are hand-curated into GO
interpretations under `curation/cgsa/` (kept separate from the native gold to
protect its human-GOA recall numbers). See `curation/cgsa/README.md`.

## Not addressed / caveats

- Running cGSA end-to-end needs an Azure OpenAI key and is not executed here.
- The published-outputs scorer uses the authors' expert category labels; it
  measures cGSA as the authors annotated it, not an independent re-judgement.
- The GO-grounding bridge only covers GO-expressible pathways; KEGG/Reactome/
  MSigDB-named ground truth stays out of GO-ID scoring by construction.
