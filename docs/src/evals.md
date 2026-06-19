# Evals

The eval framework should answer two questions:

1. Are results statistically consistent with trusted implementations on frozen inputs?
2. Is the engine faster on workloads we care about?

## Layers

Unit tests validate local pieces: bitsets, Fisher exact p-values, loaders, and annotation propagation.

CLI integration tests validate the executable surface and TSV output shape.

Eval datasets validate end-to-end workflows. These should be explicit, versioned, and reproducible.

Benchmarks validate speed on synthetic and real workloads.

## Expression20 Starter Eval

The preferred real-data smoke panel is `evals/expression20`. It names 20 human
MSigDB/MyGeneset.info expression-derived signatures selected across:

- in vivo or primary immune cell-type contrasts;
- patient or clinical disease-state contrasts;
- vaccination and infection-response contrasts;
- cytokine-stimulated primary cells;
- drug response and resistance signatures;
- cell-line stimulation/differentiation signatures.

Fetch the current GMT snapshot:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20/generated
```

Run expression-vs-expression overlap enrichment:

```bash
genesets-rs run evals/expression20/config.yaml
```

Run Expression20 against official GO:

```bash
python3 scripts/run_disease20_go_eval.py \
  --manifest evals/expression20/sets.tsv \
  --out-dir evals/expression20_vs_go/generated \
  --eval-name expression20_vs_go \
  --description "Twenty expression-derived human MSigDB/MyGeneset.info signatures enriched against official GOA human GO annotations."
```

## Disease20 Legacy Eval

Disease20 is retained as a legacy/artificial knowledge-curation contrast. It
names 20 human Disease Ontology gene sets from MyGeneset.info. The manifest is
intentionally small and hand-reviewed, but it is not representative of
differential-expression query workflows.

Fetch the current GMT snapshot:

```bash
python3 scripts/fetch_mygeneset_eval.py \
  --manifest evals/disease20/sets.tsv \
  --out-dir evals/disease20/generated
```

This writes:

- `queries.gmt`: disease gene sets in GMT format;
- `metadata.json`: source URL, fetch time, gene counts, and skipped IDs.

Generated files are not required for `cargo test`. The goal is to keep normal CI deterministic while making real-data evals one command away.

## Disease20 vs GO

The second eval uses the same 20 MyGeneset.info disease query sets against official GOA human GO annotations.

Run it:

```bash
python3 scripts/run_disease20_go_eval.py
```

This downloads official GO inputs, prepares ontology-neutral tables, runs four annotation-filter variants, and writes companion metadata:

- `all`: all evidence codes, `NOT` filtered, `contributes_to` retained;
- `no_contributes_to`: all evidence codes, `NOT` and `contributes_to` filtered;
- `iba`: IBA only, `NOT` filtered;
- `iba_iea`: IBA or IEA only, `NOT` filtered.

The prep stage writes `metadata.yaml` with source URLs, file digests, GAF header details, closure relation policy, and annotation filter counts. The run stage writes `run_metadata.yaml` with command timings and result summaries.

These filters are intentionally outside the Rust core engine. The core consumes prepared `gene_id, term_id` tables; eval helpers decide how to derive those tables from GAF.

By default, generated eval configs include `max_p_adjust: 0.05`, so result TSVs contain Bonferroni-significant rows. The companion metadata records this cutoff and the per-variant CLI runtime.

## GO Timepoints

The `go_timepoints` eval runs Disease20 against official GO archive releases:

- `2026-03-25`;
- `2021-05-01`;
- `2016-05-01`.

Run it:

```bash
python3 scripts/run_go_timepoints_eval.py
```

The summary metadata contains per-timepoint ontology sizes, GAF metadata, per-variant row counts, run timing, and pairwise diffs over `(query_id, target_id)` result keys.

## Expression500 GO All vs IBA

The IBA-impact eval compares current GOA human all-evidence annotations against
the IBA-only annotation variant over the same 500 expression-derived query sets
used by the GO impact report:

```bash
genesets-workflows go-impact evals/go_iba_impact_expression500.yaml
```

This uses the same prepared GO ontology, closure, background, and query sets on
both sides. Only the annotation table changes:

- left side: `all`, NOT-filtered GOA human, all evidence codes,
  `contributes_to` retained;
- right side: `iba`, NOT-filtered GOA human, IBA evidence only.

The local June 2026 run used the current `2026-03-25` GO/GOA prep. The
annotation variant sizes were:

| Variant | Direct gene-term pairs | Annotated genes |
| --- | ---: | ---: |
| all | 334,461 | 38,815 |
| IBA | 67,731 | 17,220 |

At Bonferroni-adjusted `p <= 0.05`, the report produced:

| Result set | Significant rows | Query sets with hits | GO terms |
| --- | ---: | ---: | ---: |
| all evidence | 26,138 | 500 | 1,125 |
| IBA only | 11,641 | 500 | 523 |

The threshold-crossing diff classified:

| Class | Rows | Query sets | GO terms |
| --- | ---: | ---: | ---: |
| lost under IBA | 17,651 | 499 | 1,014 |
| retained by IBA | 8,487 | 500 | 452 |
| gained under IBA | 3,154 | 493 | 233 |

The same run now records GO term coverage, which asks which propagated GO
targets are scorable but never become significant for this query collection:

| Annotation variant | GO terms | Scorable terms | Significant terms | Scorable but never significant | Unscorable terms |
| --- | ---: | ---: | ---: | ---: | ---: |
| all evidence | 38,560 | 21,050 | 1,125 | 19,925 | 17,510 |
| IBA only | 38,560 | 11,479 | 523 | 10,956 | 27,081 |

The paired coverage comparison classified 10,338 terms as scorable but never
significant in both variants, 9,536 as scorable but never significant on only
one side, and 1,176 as significant in at least one variant. This is useful
ontology "dark matter" for a fixed query collection: terms can be fully
annotated and still never show up in enrichment output.

This supports the expected contraction from all evidence to IBA-only, but it
also shows why post-processing matters. Because IBA changes target sizes and
the correction universe, some very broad terms become significant only in the
IBA run. With a simple specificity guard of target size `<= 1,000`, the diff
becomes 3,699 lost, 1,401 retained, and 350 gained rows.

Interpretation: IBA is a useful vetted subset, but not an automatic
"less-biased answer set." It is also a lower-coverage subset, and it can shift
enrichment toward terms whose IBA target sizes shrink dramatically. Reports
should therefore present all-vs-IBA threshold crossings together with target
size deltas and a configurable slim/antislim or size-based term filter.

Qualitative spot checks support that caution:

| Query set | Observation |
| --- | --- |
| `GSE29618_PDC_VS_MDC_DAY7_FLU_VACCINE_DN` | All-evidence losses include `secretory granule` and `regulation of cytokine production`, with overlaps dropping from 39 to 1 and 32 to 6 respectively under IBA. These terms make sense for dendritic-cell/vaccine biology; IBA is mainly losing immune-cell annotation coverage. |
| `GSE30962_PRIMARY_VS_SECONDARY_ACUTE_LCMV_INF_CD8_TCELL_UP` | All-evidence losses include cell-cycle regulation terms supported by genes such as `AURKB`, `BIRC5`, `BUB1`, `CDC6`, and `RRM2`. The biology fits proliferating/secondary CD8 T-cell response. IBA retains too little of that term coverage. IBA-only gains such as meiotic recombination are driven by the same DNA repair/recombination genes and are semantically less appropriate for CD8 T cells. |
| `GSE31082_DN_VS_DP_THYMOCYTE_UP` | All-evidence mitochondrial-envelope losses are supported by many mitochondrial/ribosomal genes. IBA-only gains emphasize RNA metabolic and ribonucleoprotein terms because those target sets shrink sharply; this is plausible for thymocyte transition but is a reframing rather than a new signal. |
| `GSE11924_TFH_VS_TH17_CD4_TCELL_DN` | IBA-only gains are mostly ribonucleoprotein/RNA-processing terms. They fit a generic translation/RNA-processing axis, not a helper-T-cell-specific interpretation. |
| `GSE13485_CTRL_VS_DAY21_YF17D_VACCINE_PBMC_DN` | IBA-only synapse/cell-junction gains are supported by genes such as cholinergic/glutamate receptor genes and synaptic adhesion genes. This may be real in the input gene set, but it is not an obvious vaccine/PBMC biology headline and should be treated cautiously. |
| `GSE28726_NAIVE_VS_ACTIVATED_CD4_TCELL_DN` | All-evidence losses include cell-cycle regulation terms supported by canonical proliferation genes such as `AURKA`, `AURKB`, `BRCA1`, `CDC20`, `CDC6`, `CENPF`, `MKI67`, and `UBE2C`. These losses make biological sense for an activation contrast; IBA is filtering away coverage rather than clearly removing bias. |

The working interpretation is therefore: IBA is more curated and useful as a
contrast set, but many "lost" calls are plausible biology that lacks IBA
coverage. Many "gained" calls are not discoveries; they are threshold crossings
caused by smaller IBA target sets and a changed correction universe. The report
should expose target-size deltas and representative overlap genes before
promoting a gained/lost term as biologically meaningful.

The losses are not all equivalent. Useful report examples should separate
desirable losses from concerning losses:

Good losses, where IBA removes broad or weakly interpretable all-evidence
calls:

| Lost term | Why this is probably good |
| --- | --- |
| `protein binding` | Lost in 449 query sets with an all-evidence target size of 15,319 genes. It is usually too broad to be useful as a headline enrichment term. |
| `binding` | Lost in 238 query sets with an all-evidence target size of 19,997 genes. This is a classic uninformative molecular-function call. |
| `cytosol` | Lost in 377 query sets. Often useful as context, but as a recurring top-level cellular-component hit it can dominate reports without adding much interpretation. |
| `membrane` / `endomembrane system` / `vesicle` | These are often real but broad compartment signals. Dropping some of them can improve report readability, especially when more specific descendants remain. |
| `positive regulation of biological process` | A very broad biological-process regulation call; loss is usually desirable unless a more specific regulation term is also lost. |

Bad losses, where IBA removes biologically plausible calls for the query set:

| Query set | Lost term | Why this is concerning |
| --- | --- | --- |
| `GSE29618_PDC_VS_MDC_DAY7_FLU_VACCINE_DN` | `secretory granule` | Strong all-evidence overlap drops from 39 query genes to 1 IBA-covered overlap. This fits dendritic/myeloid biology and looks like coverage loss. |
| `GSE29618_PDC_VS_MDC_DAY7_FLU_VACCINE_DN` | `regulation of cytokine production` | Overlap drops from 32 to 6. This is directly relevant to vaccine/dendritic-cell biology. |
| `GSE30962_PRIMARY_VS_SECONDARY_ACUTE_LCMV_INF_CD8_TCELL_UP` | `positive regulation of cell cycle process` | Supported by proliferation genes such as `AURKB`, `BIRC5`, `BUB1`, `CDC6`, and `RRM2`. This is credible CD8 response biology. |
| `GSE28726_NAIVE_VS_ACTIVATED_CD4_TCELL_DN` | `regulation of cell cycle process` | Supported by `AURKA`, `AURKB`, `BRCA1`, `CDC20`, `CDC6`, `CENPF`, `MKI67`, and `UBE2C`. This makes sense for activation state. |
| `GSE31082_DN_VS_DP_THYMOCYTE_UP` | `mitochondrial envelope` | Overlap drops from 40 to 7. The query contains many mitochondrial/ribosomal genes, so this looks like a real thymocyte-transition signal. |

The report should therefore label losses as "possibly useful pruning" only
when they are broad, recurrent, and replaceable by more specific retained
terms. Losses that match the experimental contrast and have coherent overlap
genes should be highlighted as IBA coverage gaps.

Some bad losses have partial compensation, but the compensation is usually from
retained IBA-significant neighbors rather than from newly gained IBA-only
terms:

| Query set | Lost term | Compensatory signal |
| --- | --- | --- |
| `GSE29618_PDC_VS_MDC_DAY7_FLU_VACCINE_DN` | `regulation of cytokine production`, `secretory granule` | IBA retains immune terms such as `immune system process`, `regulation of immune system process`, `positive regulation of immune system process`, and lipid-antigen presentation terms. IBA-only gains include lipid-antigen binding, but the only ontology ancestors gained for the lost terms are very broad terms such as `biological_process` and `intracellular organelle`, so they are not good semantic replacements. |
| `GSE30962_PRIMARY_VS_SECONDARY_ACUTE_LCMV_INF_CD8_TCELL_UP` | `positive regulation of cell cycle process` | Strong compensation exists through retained IBA terms: `cell cycle`, `cell cycle process`, `mitotic cell cycle`, `chromosome segregation`, `nuclear division`, and `DNA replication`. The IBA-only gains are more generic or less cell-type-appropriate, such as nucleic-acid metabolism and meiotic recombination. |
| `GSE28726_NAIVE_VS_ACTIVATED_CD4_TCELL_DN` | `regulation of cell cycle process` | IBA still retains `cell cycle`, `cell cycle process`, `mitotic cell cycle`, `DNA replication`, and related chromosome terms. The gained terms again skew toward nucleic-acid metabolism and meiotic/catalytic terms, which are weaker replacements. |
| `GSE31082_DN_VS_DP_THYMOCYTE_UP` | `mitochondrial envelope` | IBA retains nearby mitochondrial terms such as `mitochondrion`, `mitochondrial matrix`, `mitochondrial protein-containing complex`, and `mitochondrial ribosome`. IBA-only gains emphasize RNA metabolism, ribosome, and translation, which may be a real thymocyte-transition axis but do not fully replace the envelope/membrane interpretation. |

These spot checks are still too immune-heavy. That is partly a fixture problem:
the first 500 MyGeneset/MSigDB `GSE*` hits are dominated by C7-style
immunologic signatures. The IBA report should therefore avoid presenting those
examples as representative of all expression data. Less immune-centric checks
from the same run are more mixed:

| Query set | Observation |
| --- | --- |
| `GSE11367_CTRL_VS_IL17_TREATED_SMOOTH_MUSCLE_CELL_DN` | The lost calls are broad terms such as `response to stimulus`, `extracellular region`, `vesicle`, `membrane`, and `protein binding`. The biology is plausible for IL17-treated smooth muscle, but the terms are too general to make this a strong "bad loss" example. |
| `LI_STAD_HAZARD_RATIO_HIGH` | No significant all-evidence calls are lost. IBA gains are broad organelle/root-level terms. This is useful as a negative example: not every threshold crossing deserves interpretation. |
| `GSE12003_MIR223_KO_VS_WT_BM_PROGENITOR_4D_CULTURE_UP` | Non-immune-specific losses include `nuclear speck` and `chromatin binding`. `nuclear speck` drops from 19 overlapping query genes in all evidence to 1 under IBA; `chromatin binding` drops from 20 to 6. These look like real coverage losses around nuclear/chromatin biology. |
| `GSE12003_4D_VS_8D_CULTURE_BM_PROGENITOR_UP` | `ribosome biogenesis` is lost, dropping from 14 overlapping genes to 9 under IBA and falling below the cutoff. IBA gains `nucleolus`, with 12 overlapping IBA genes, so there is partial compensation in a neighboring cellular-component term. |
| `GSE29164_DAY3_VS_DAY7_UNTREATED_MELANOMA_DN` | All-evidence `endoplasmic reticulum` and `endoplasmic reticulum membrane` calls are lost, dropping from 46 to 17 and 30 to 8 overlapping genes respectively. IBA still retains `endomembrane system`, so the broad compartment signal remains but the ER-specific interpretation disappears. |
| `GSE29164_DAY3_VS_DAY7_UNTREATED_MELANOMA_UP` | The specific lost calls are lymphocyte/T-cell activation terms. This is probably not melanoma-cell-intrinsic biology; it is a reminder that tumor-labeled signatures can carry immune microenvironment signal. |

This suggests a useful display rule: for each lost term, show whether the same
query still has retained or gained terms in the same ontology neighborhood. A
loss with strong neighboring retained terms is less damaging than a loss where
the entire biological neighborhood disappears.

For the next expression-scale fixture, sample by source family or keywords
rather than taking the first `GSE*` results. The eval set should deliberately
include immune, cancer, stromal, epithelial, neuronal, developmental,
cell-cycle, drug-response, and stress-response signatures so qualitative
examples are not all drawn from immune contrasts.

The stratified 5,000-set plan lives at
`evals/expression_like/msigdb_diverse_5k.yaml` and can be fetched with:

```bash
genesets-workflows fetch-mygeneset-stratified \
  evals/expression_like/msigdb_diverse_5k.yaml
```

The same source can drive temporal GO and all-vs-IBA reports:

```bash
genesets-workflows go-impact evals/go_impact_5y_expression5000_diverse.yaml
genesets-workflows go-impact evals/go_iba_impact_expression5000_diverse.yaml
```

Benchmark composition matters. The refreshed source reserves modern MSigDB
families rather than relying on first-N `GSE*` or old pathway hits:

| Source family | Sets |
| --- | ---: |
| C8-like single-cell marker signatures | 496 |
| C9 DepMap/CCLE perturbation signatures | 62 |
| C4/3CA cancer metaprograms | 148 |
| GSE expression signatures | 1,212 |
| GO-derived controls | 687 |
| HPO phenotype-derived sets | 313 |
| Pathway-derived sets | 362 |
| Hallmark sets | 12 |
| Other curated MSigDB sets | 1,708 |

This is a benchmark-suite rule, not just a convenience for one notebook:
large gene-set benchmarks should report source-family composition and reserve
space for newer collections such as MSigDB C8, C9, and 3CA.

A local synthetic scale probe used the current 500-set fixture repeated with
unique query IDs. This is not a biological eval, but it answers the engineering
question about 10x query-set scale:

| Workload | Runtime |
| --- | ---: |
| 500 query sets, median 209 genes, null output | 4.3s |
| 500 query sets, median 1,772 genes, null output | 4.9s |
| 5,000 query sets, null output | 15.4s |
| 5,000 query sets, Parquet output | 15.8s |
| 5,000 query sets, IBA-only Parquet output | 14.9s |
| 5,000 all-vs-IBA Parquet compare | 1.2s |

The full 5,000-set synthetic all-vs-IBA path was about 31s on this laptop. The
all-evidence side wrote 208,220 significant rows; the IBA side wrote 99,290;
the diff wrote 237,560 rows. The row count is not exactly 10x the 500-set run
because Bonferroni becomes stricter as query count increases.

The actual stratified 5,000-set all-vs-IBA report completed in 52.6s after
refreshing the query source to include explicit C8/C9/3CA coverage:

| Step | Runtime | Rows |
| --- | ---: | ---: |
| all-evidence matrix | 23.9s | 473,603 |
| IBA-only matrix | 24.3s | 163,866 |
| compare | 3.2s | 501,223 |

The diff classes were 337,357 lost, 136,246 shared, and 27,620 gained
significant pairs. This is a better scale signal than the repeated synthetic
run because the query collection is actually diverse and has a broader gene
background.

The same stratified source against the five-year temporal GO/GOA diff completed
in 47.3s:

| Step | Runtime | Rows |
| --- | ---: | ---: |
| 2021-05-01 matrix | 17.7s | 289,007 |
| 2026-03-25 matrix | 24.7s | 473,603 |
| compare | 3.7s | 546,972 |

The temporal diff classes were 73,369 lost, 215,638 shared, and 257,965 gained
significant pairs. This is the current best estimate for a coarse two-snapshot
GO impact report over 5,000 diversified expression-like query sets.

## Expression-Like Query Sets

Disease Ontology gene sets are convenient but artificial. For scale and realism, use MyGeneset/MSigDB signatures:

```bash
genesets-workflows fetch-mygeneset \
  --query 'GSE*' \
  --source-filter msigdb \
  --limit 2000 \
  --out-dir evals/expression_like/generated/msigdb_gse_2k
```

MSigDB GSE, `*_UP`, `*_DN`, and `vs` signatures are closer to differential-expression outputs than disease ontology gene sets. These should become the main large-query workloads for timing and GO-version diffing.

## Expression500 vs Reactome Flat

Reactome can be tested without any new Rust code because the core engine
already supports flat GMT target libraries. The current flat Reactome eval uses
the same 500 expression-derived query sets as the GO impact report and enriches
them against the official Reactome pathway GMT.

Run it:

```bash
genesets-workflows reactome-flat
```

This performs two steps:

1. `genesets-workflows prepare-reactome-flat` downloads
   `ReactomePathways.gmt.zip` from the
   official Reactome download directory and normalizes it from `name, id,
   genes...` to `id, name, genes...`.
2. `genesets-rs run evals/reactome_flat/config.yaml` enriches the expression
   query sets against that flat pathway library.

The local May 2026 run prepared 2,557 Reactome pathways over 11,963 gene
symbols and produced 509 Bonferroni-significant rows across 85 query sets. The
enrichment step took about 0.22 seconds on this laptop.

This is intentionally the flat baseline. A hierarchy-aware Reactome workflow
should prepare the standard ontology-style trio:

- `terms.tsv` from `ReactomePathways.txt`;
- `closure.tsv` from `ReactomePathwaysRelation.txt`;
- `gene_terms.tsv` from a lowest-level mapping such as `NCBI2Reactome.txt`,
  with a gene ID strategy matched to the query sets.

## Planned Comparison Evals

The most useful comparison suite will freeze:

- ontology release;
- association release;
- closure relation policy;
- gene universe;
- query gene sets;
- exact correction method;
- expected significant terms from each comparator.

Comparator adapters should run each external tool in containers where possible and emit normalized TSV for comparison.
