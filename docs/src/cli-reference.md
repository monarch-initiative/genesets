# CLI Reference

The CLI has four subcommands:

- `enrich`: one query set against many target sets;
- `matrix`: many query sets against many target sets;
- `run`: YAML-configured `enrich` or `matrix`.
- `compare`: threshold-crossing diff between two result tables.

## `enrich`

```bash
genesets-rs enrich \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --sample sample.txt \
  --background background.txt \
  --output results.tsv
```

Use `--target-sets` instead of `--annotations` for flat libraries:

```bash
genesets-rs enrich \
  --target-sets library.gmt \
  --target-format gmt \
  --sample sample.txt
```

Useful options:

- `--sample-format`: `auto`, `list`, `pairwise`, `gene-term`, `gmt`, `gmx`, `gmx-desc`;
- `--sample-set`: select one set from a multi-set sample file;
- `--min-overlap`: suppress rows below an overlap count;
- `--max-p-value`: suppress rows above a raw p-value cutoff;
- `--max-p-adjust`: suppress rows above an adjusted p-value cutoff;
- `--correction`: `bonferroni` or `none`;
- `--output-format`: `tsv`, `parquet`, or `null`; `parquet` requires `--output`,
  and `null` is useful for compute-only profiling;
- `--overlap-genes`: include overlapping gene IDs and names;
- `--threads`: set Rayon worker count.

## `matrix`

```bash
genesets-rs matrix \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --queries queries.gmx \
  --query-format gmx \
  --background background.txt
```

Use targets as queries for term-vs-term runs:

```bash
genesets-rs matrix \
  --annotations gene_terms.tsv \
  --closure closure.tsv \
  --queries-from-targets
```

## `run`

```yaml
mode: enrich
ontology:
  terms: terms.tsv
  closure: closure.tsv
  annotations: gene_terms.tsv
input:
  sample: sample.txt
  sample_format: list
  sample_name: sample
background:
  file: background.txt
overlap_genes: true
max_p_adjust: 0.05
output_format: tsv
```

Run it:

```bash
genesets-rs run examples/enrich.yaml
```

Relative paths in YAML are resolved from the config file's directory. Relative
paths passed directly as CLI arguments are resolved from the current working
directory.

For mass evals, write Parquet and inspect it with DuckDB:

```bash
genesets-rs matrix ... --output-format parquet --output results.parquet
duckdb -c "SELECT * FROM 'results.parquet' WHERE p_adjust_bonferroni <= 0.05"
```

## `mgsa`

Bayesian model-based gene set analysis. It takes the same target and sample
options as `enrich`. Every set in a multi-set sample file is analyzed
independently unless `--sample-set` picks one.

```bash
genesets-rs mgsa \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --sample study_genes.txt \
  --background population_genes.txt \
  --min-posterior 0.5 \
  --summary-output mgsa-summary.yaml
```

Output columns: `query_id`, `query_name`, `target_id`, `target_name`,
`posterior`, `overlap`, `query_size`, `target_size`, `background_size`,
`p_value` (term-for-term, for comparison), `equivalent_ids`, and with
`--overlap-genes`, `overlap_genes` and `overlap_gene_names`. Rows are sorted by
posterior, highest first.

| Option | Default | Meaning |
| --- | --- | --- |
| `--samples-per-term` | 100 | recorded samples per modeled term |
| `--burn-per-term` | 10 | burn-in samples per modeled term |
| `--term-prob`, `--term-count` | 0.5, 0 | Beta prior on a term being active |
| `--false-pos-prob`, `--false-pos-count` | 0.5, 0 | Beta prior on the false-positive rate |
| `--false-neg-prob`, `--false-neg-count` | 0.5, 0 | Beta prior on the false-negative rate |
| `--flip-rate`, `--step-rate`, `--jump-rate` | 1, 1, 1 | relative move rates; `step` needs `--closure` |
| `--seed` | 123456789 | random seed |
| `--min-posterior` | 0 | drop rows below this posterior |

Only TSV (or `null`) output is supported. See [Statistics](statistics.md) for
the model.

## `gsea`

Preranked gene set enrichment analysis. It takes the same target options as
`enrich`, but a scored gene list instead of a sample set and background: the
ranking is the universe, so every ranked gene counts whether or not it is
annotated, and target sets are restricted to ranked genes.

```bash
genesets-rs gsea \
  --annotations gene_terms.tsv \
  --terms terms.tsv \
  --closure closure.tsv \
  --ranks contrast.rnk \
  --max-p-adjust 0.05 \
  --leading-edge-genes
```

`--ranks` accepts two columns (`gene`, `score`), as in a GSEA `.rnk` file, or
three (`ranking_id`, `gene`, `score`) to run several rankings at once; each is
analyzed and adjusted independently. A first row with a non-numeric score is
treated as a header. Duplicate genes and non-finite scores are errors.

Output columns: `query_id` (the ranking), `query_name`, `target_id`,
`target_name`, `es`, `nes`, `p_value`, `p_adjust` (Benjamini-Hochberg within a
ranking), `set_size` (members in the ranking), `target_size`, `ranked_size`,
`leading_edge_size`, `peak_rank` (1-based rank where the running sum peaks), and
with `--leading-edge-genes`, `leading_edge_genes` and
`leading_edge_gene_names`. Rows are sorted by p-value, then by |NES|. The
`p_adjust` column means `compare` can diff two GSEA runs directly.

| Option | Default | Meaning |
| --- | --- | --- |
| `--ranks-name` | file stem | ranking id for a two-column file |
| `--ranks-set` | all | analyze one ranking from a three-column file |
| `--weight` | 1 | exponent on \|score\|; 0 is the unweighted KS statistic |
| `--permutations` | 1000 | random gene sets per distinct set size |
| `--min-size`, `--max-size` | 15, 500 | ranked-member size limits for tested targets |
| `--seed` | 123456789 | random seed |
| `--max-p-value`, `--max-p-adjust` | none | drop rows above these cutoffs |

Only TSV (or `null`) output is supported. See [Statistics](statistics.md#preranked-gsea)
for the method.

## `compare`

Compare two enrichment result tables by `(query_id, target_id)` and classify
adjusted p-value threshold crossings:

```bash
genesets-rs compare \
  --left go-2021.parquet \
  --right go-2026.parquet \
  --p-adjust-cutoff 0.05 \
  --output-format parquet \
  --output go-2021-vs-2026.diff.parquet \
  --metadata-output go-2021-vs-2026.diff.yaml
```

Input formats are inferred from `.tsv`, `.txt`, `.parquet`, or `.pq`, or can be
set explicitly with `--left-format` and `--right-format`. TSV output goes to
stdout by default. Parquet output requires `--output`.

The default output includes:

- `lost_significant`: significant on the left, not significant on the right;
- `gained_significant`: not significant on the left, significant on the right;
- `shared_significant`: significant on both sides.

Use `--crossings-only` to emit only gained/lost rows.

## Workflow CLI

The separate `genesets-workflows` command is the convenience layer for
configured source prep and reports. It calls the Rust CLI for batch compute and
then writes Parquet plus metadata:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
genesets-workflows go-impact evals/go_impact_5y_expression500.yaml
genesets-workflows reactome-flat
```

`go-impact` YAML can apply `query_sets.include_regex` and
`query_sets.exclude_regex` before selecting query sets from a cached snapshot.
The regexes match query id, name, description, source class, stratum, and search
query metadata. For example, exclude GO-derived MSigDB query sets from a GO
target analysis:

```yaml
query_sets:
  source_dir: expression_like/generated/msigdb_diverse_5k
  limit: 4313
  exclude_regex:
    - "^(GOBP|GOCC|GOMF)_"
```

For interactive review of an existing report bundle, use the optional web
explorer dependencies:

```bash
uv run --project python/genesets-workflows --extra explorer \
  genesets-workflows explore notebooks/generated/go_iba_impact_expression5000_diverse
```

Use `genesets-rs` for normalized single jobs. Use `genesets-workflows` when the
task needs downloads, evidence filters, release metadata, multiple Rust runs,
DuckDB summaries, notebook/report artifacts, or browser-based result triage.
