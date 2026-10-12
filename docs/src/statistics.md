# Statistics

`enrich` and `matrix` perform standard one-sided over-representation analysis.
`mgsa` fits a Bayesian model over all terms jointly; see
[Model-Based Gene Set Analysis](#model-based-gene-set-analysis).

For each query and target term, the engine computes:

| Symbol | Meaning |
| --- | --- |
| `N` | background size |
| `K` | target size in background |
| `n` | query size in background |
| `k` | query-target overlap in background |

The p-value is the right-tail hypergeometric probability:

```text
P(X >= k), X ~ Hypergeometric(N, K, n)
```

This is equivalent to one-sided Fisher exact enrichment on the 2 x 2 table:

| | In target | Not in target |
| --- | ---: | ---: |
| In query | `k` | `n - k` |
| Not in query | `K - k` | `N - K - n + k` |

## Multiple Testing

The current correction choices are:

- `bonferroni`: `min(p * number_of_tests, 1)`;
- `none`: report raw p-values as adjusted p-values.

For matrix runs, the Bonferroni denominator is the number of non-empty query by non-empty target tests.

## Model-Based Gene Set Analysis

`genesets-rs mgsa` implements the Bayesian model of Bauer, Gagneur & Robinson
(2010), [GOing Bayesian: model-based gene set analysis of genome-scale
data](https://doi.org/10.1093/nar/gkq045), which Ontologizer calls MGSA. The
implementation follows the collapsed form in Ian Holmes and Chris Mungall's
[wtfgenes](https://github.com/evoldoers/wtfgenes).

Term-for-term tests score each term on its own, so a strongly enriched term
drags its ancestors and overlapping neighbors up with it. MGSA asks a
different question: which small set of *active* terms best explains the whole
query?

- Each term is active or inactive. A gene annotated to at least one active term
  is "on"; other genes are "off".
- An off gene appears in the query with false-positive rate `alpha`. An on gene
  is missing from the query with false-negative rate `beta`. Each term is active
  a priori with rate `p`.
- `p`, `alpha`, and `beta` each have a Beta prior and are integrated out
  analytically, so the likelihood depends only on six counts: active and
  inactive terms, off genes in and out of the query, and on genes in and out of
  the query.
- A Metropolis-Hastings sampler walks over term states. `flip` toggles one term,
  `step` swaps an active term for a parent, child, or sibling, and `jump` swaps
  an active term for any modeled term. `step` needs `--closure`; it is skipped
  for flat target sets.

The reported `posterior` is the fraction of post-burn-in samples in which the
term was active. Like wtfgenes, the run length scales with the number of modeled
terms: `--samples-per-term` (default 100) and `--burn-per-term` (default 10).

Which terms are modeled:

- only terms with at least one query gene in the background;
- terms with identical genes in the background form one variable. The term with
  the most ancestors speaks for the group, and the others are listed in
  `equivalent_ids`.

Priors are set by mode and pseudocount, as in wtfgenes: `--term-prob` and
`--term-count`, `--false-pos-prob` and `--false-pos-count`, `--false-neg-prob`
and `--false-neg-count`. The defaults (mode 0.5, zero pseudocounts) are uniform
priors. `--summary-output` writes per-query diagnostics as YAML: modeled terms,
sample counts, acceptance rate, mean number of active terms, and the posterior
means of `p`, `alpha`, and `beta`.

The sampler is seeded (`--seed`), so runs are reproducible. Each query in a
multi-set sample file is analyzed independently and in parallel. Posteriors
are Monte Carlo estimates. If a term's posterior changes noticeably between
seeds, raise `--samples-per-term`.

Run times and a comparison with wtfgenes are in
[Performance Model](performance.md#mgsa).

Differences from Ontologizer's MGSA: Ontologizer samples `alpha`, `beta`, and
`p` over a discrete grid alongside the term states, while this implementation
integrates them out against Beta priors. Expect similar rankings rather than
identical posteriors. On the `examples/ontologizer-ui` study and population
lists with current human GOA, all nine terms that the Ontologizer UI reports as
high-scoring rank in the top 12 here.

## Future Statistics

The code is structured so additional test families can be added without changing the data model:

- Benjamini-Hochberg FDR;
- ranked-list methods;
- topology-aware parent-child methods;
- elim and weight-style algorithms.
