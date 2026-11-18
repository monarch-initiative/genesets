# Performance Model

`genesets-rs` trades initialization work and memory for fast repeated calculations.

## Bitset Representation

Each query and target is represented as a dense `Vec<u64>` over a shared gene universe. After loading:

- set size is a popcount;
- overlap is a word-wise AND plus popcount;
- ontology propagation is already baked into target bitsets;
- matrix mode parallelizes over query sets.

This is a good fit when the universe is moderate, memory is available, and the same target library is reused across many query sets.

## Expected Bottlenecks

There are three separate costs to keep apart:

- preparing inputs: parsing term, closure, annotation, background, and query files, then building dense bitsets;
- scoring: intersecting query and target bitsets, running Fisher exact tests, and applying cutoffs;
- writing results: resolving string labels and serializing rows.

The CLI has `--output-format null` for compute-only timing. Use it to separate scoring from serialization:

```bash
genesets-rs matrix ... --output-format null
```

On the current 200 MSigDB GSE-style query fixture against GOA human and the current GO closure:

- one query, significant-only, no output serialization: about 1.4s with one Rayon thread;
- 200 queries, significant-only, no output serialization: about 5.9s with one Rayon thread, or about 2.2s on a 10-core laptop;
- 200 queries, no p-value cutoff, no output serialization: about 2.25s wall time and about 540 MB peak RSS;
- 200 queries, no p-value cutoff, TSV output: about 3.6s wall time for about 778k result rows and a 140 MB TSV.

For significant-only evals, TSV is not the main bottleneck. Rebuilding the ontology target index on every CLI process dominates wall time. For unfiltered or weakly filtered mass comparisons, row materialization and TSV serialization become large enough to matter.

TSV remains useful for individual analyses, but multi-timepoint and many-query evals should move to a columnar or embedded analytical format. The output module is intentionally isolated from the scoring layer so Parquet and DuckDB sinks can be added without changing enrichment statistics.

## Current Diagnosis

The most important bottlenecks are:

- repeated text parsing and ontology target bitset construction for every run;
- duplicated string-heavy closure and annotation tables during setup;
- materializing every retained result row in memory before writing;
- sorting all retained rows before output;
- dense bitset scans over every query-target pair, even for very sparse terms;
- TSV formatting when result sets are large or unfiltered.

The next performance milestones are:

- serialized prepared indexes for term bitsets, term metadata, and gene metadata;
- streaming or batched result sinks so large runs do not materialize every row;
- Parquet output for columnar scans and archival eval artifacts;
- DuckDB output for direct diff queries over many timepoints;
- top-k and p-value-threshold output modes;
- mmap-friendly bitset storage;
- sparse or Roaring-style bitsets benchmarked against dense bitsets on real GO fixtures;
- criterion benchmarks on synthetic and real ontology fixtures;
- comparison runs against GOATOOLS, g:Profiler, topGO, and clusterProfiler on the same frozen inputs.

## MGSA

`genesets-rs mgsa` is a Markov chain Monte Carlo (MCMC) sampler, so it costs
more than Fisher scoring. The run length is
`(samples_per_term + burn_per_term) x modeled_terms` steps (110 per term by
default). Each step touches only the genes of the one or two terms it toggles
and updates six integer counts. The sampler reaches roughly 200k-500k steps per
second on one core; the rate depends on term sizes and the acceptance rate.

Timings from a 4-core cloud container with current GO (is_a plus part_of) and
human GOA. Loading GO and GOA takes about 1.7s of each run and is the same cost
as for `enrich`:

| Job | Modeled terms | Wall time |
| --- | ---: | ---: |
| `enrich` on the `examples/ontologizer-ui` lists (load-only baseline) | n/a | 1.6-1.7s |
| `mgsa`, same lists (969 study genes, 19,687 population genes) | 8,740 | 3.4-4.1s |
| `mgsa`, same study genes, background = every GOA gene, 1 thread | 8,759 | 6.1s |

Sampling cost per query grows with query size, because larger queries touch
more terms. These are random population genes, one thread, load time excluded:

| Query size | Modeled terms | Sampling per query |
| ---: | ---: | ---: |
| 50 | ~1,500 | ~0.6s |
| 200 | ~3,500 | ~0.85s |
| 1,000 | ~7,500 | ~1.4s |
| 3,000 | ~11,400 | ~2.0s |

Queries in a multi-set sample file run in parallel: 8 queries of 3,000 genes
took 17.4s on one thread and 6.2s on four.

For reference, the original JavaScript wtfgenes (`bin/wtfgenes.js`, defaults)
ran the 6.1s job above at about 500 steps per second. It reached 32% of its
958k steps after a 10-minute timeout, and its progress meter projected about
30 minutes in total. The gap comes mostly from data structures rather than the
language. wtfgenes keeps counts in string-keyed objects and allocates new ones
on every proposal, while genesets-rs updates dense arrays in place.

To trade accuracy for speed, lower `--samples-per-term`; run time scales
roughly linearly with it. Use `--summary-output` to see the modeled-term count,
the acceptance rate, and the number of samples behind each posterior.

## Result Storage

For the mass-diff use case, the preferred storage shape is numeric and columnar:

- dimension tables for runs, queries, targets, terms, genes, ontology releases, and annotation variants;
- fact rows keyed by `(run_id, query_index, target_index)`;
- numeric columns for overlap, sizes, p-value, adjusted p-value, and significance flags;
- optional overlap-gene payloads in a separate table or sidecar file.

Parquet is a good artifact format for immutable eval output and batch scans. DuckDB is a good working format for interactive diff queries, threshold-crossing queries, joins against ontology metadata, and summaries across many releases. Both should sit behind the output boundary rather than leak into the enrichment model.

## Benchmark Harness

The repository includes a Criterion benchmark:

```bash
cargo bench
```

The benchmark is synthetic by design. It measures the core scoring path without network access, identifier mapping, or ontology parsing.
