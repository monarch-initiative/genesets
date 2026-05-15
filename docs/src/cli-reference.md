# CLI Reference

The CLI has three subcommands:

- `enrich`: one query set against many target sets;
- `matrix`: many query sets against many target sets;
- `run`: YAML-configured `enrich` or `matrix`.

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
  terms: examples/terms.tsv
  closure: examples/closure.tsv
  annotations: examples/gene_terms.tsv
input:
  sample: examples/sample.txt
  sample_format: list
  sample_name: sample
background:
  file: examples/background.txt
overlap_genes: true
max_p_adjust: 0.05
output_format: tsv
```

Run it:

```bash
genesets-rs run examples/enrich.yaml
```

For mass evals, write Parquet and inspect it with DuckDB:

```bash
genesets-rs matrix ... --output-format parquet --output results.parquet
duckdb -c "SELECT * FROM 'results.parquet' WHERE p_adjust_bonferroni <= 0.05"
```
