# Expression-Like Query Evals

Disease Ontology sets are useful smoke tests, but they are not realistic query sets for expression-analysis workflows.

For the small curated panel, use `evals/expression20`. This directory is for
larger MyGeneset/MSigDB query snapshots that look more like
differential-expression signatures.

Examples:

```bash
genesets-workflows fetch-mygeneset \
  --query 'GSE*' \
  --source-filter msigdb \
  --limit 2000 \
  --out-dir evals/expression_like/generated/msigdb_gse_2k
```

```bash
genesets-workflows fetch-mygeneset \
  --query 'vs' \
  --source-filter msigdb \
  --limit 5000 \
  --out-dir evals/expression_like/generated/msigdb_vs_5k
```

These query sets can then be used with GO prepared configs by replacing `input.queries` with the generated `queries.gmt`.

The first-N `GSE*` snapshots are useful for quick scale checks but are biased
toward immunologic signatures. For biological spot checks and report examples,
prefer the stratified plan:

```bash
genesets-workflows fetch-mygeneset-stratified \
  evals/expression_like/msigdb_diverse_5k.yaml
```

This writes `evals/expression_like/generated/msigdb_diverse_5k/queries.gmt`
with per-stratum provenance in `metadata.json`. The current quota plan targets
5,000 MSigDB sets and reserves explicit space for newer MSigDB families
described by the official MSigDB collection page, including
[C8 cell-type signatures, C9 perturbation signatures, and C4/3CA cancer
metaprograms](https://www.gsea-msigdb.org/gsea/msigdb/collections.jsp), before
filling broad keyword strata:

| Reserved/source family | Sets in refreshed source |
| --- | ---: |
| C8-like single-cell marker signatures | 496 |
| C9 DepMap/CCLE perturbation signatures | 62 |
| C4/3CA cancer metaprograms | 148 |
| GSE expression signatures | 1,212 |
| GO-derived MSigDB controls | 687 |
| HPO phenotype-derived sets | 313 |
| Pathway-derived sets | 362 |
| Hallmark sets | 12 |
| Other curated MSigDB sets | 1,708 |

The strata are intentionally non-perfect keyword buckets; the goal is fixture
diversity and auditable provenance, not ontology-grade sample classification.

The intended scale tests are:

- 2k signatures x GO, single timepoint;
- 10k signatures x GO, single timepoint;
- 2k signatures x 3 GO timepoints x 4 GAF variants;
- larger runs once output filtering and diff summaries are stable.

Local June 2026 synthetic scale probe against current GOA human/GO:

| Query workload | Output | Runtime |
| --- | --- | ---: |
| 500 existing expression-like sets, median 209 genes | null | 4.3s |
| 500 synthetic sets, median 1,772 genes | null | 4.9s |
| 5,000 repeated expression-like sets | null | 15.4s |
| 5,000 repeated expression-like sets | Parquet | 15.8s |
| 5,000 repeated expression-like sets, IBA only | Parquet | 14.9s |
| 5,000 all-vs-IBA compare | Parquet | 1.2s |

Interpretation: 10x more query sets scales comfortably on this laptop for the
current filtered-result workflow. Larger per-query gene counts are not the main
bottleneck with the current dense-bitset representation. The more important
future scaling risk is the number of query-set by target-term tests retained
for output and downstream diff summaries.

The actual stratified 5,000-set all-vs-IBA run completed in 52.6s after
refreshing the source to reserve C8/C9/3CA coverage:

| Step | Runtime | Significant rows |
| --- | ---: | ---: |
| all-evidence matrix | 23.9s | 473,603 |
| IBA-only matrix | 24.3s | 163,866 |
| compare | 3.2s | 501,223 diff rows |

The generated source had 5,000 sets, median 117 genes per set, 742,649 listed
gene memberships, and 24,849 unique background genes.

The same source for the 2021-05-01 vs 2026-03-25 temporal GO report completed
in 47.3s:

| Step | Runtime | Significant rows |
| --- | ---: | ---: |
| 2021-05-01 matrix | 17.7s | 289,007 |
| 2026-03-25 matrix | 24.7s | 473,603 |
| compare | 3.7s | 546,972 diff rows |
