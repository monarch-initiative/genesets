# Statistics

`enrich` and `matrix` perform standard one-sided over-representation analysis.
`mgsa` fits a Bayesian model over all terms jointly; see
[Model-Based Gene Set Analysis](#model-based-gene-set-analysis). `gsea` scores
terms against a ranked, scored gene list instead of a set; see
[Preranked GSEA](#preranked-gsea).

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

From Python, `genesets.mgsa(...)` in `python/genesets-py` runs the same
sampler in-process; see that package's README.

Run times and a comparison with wtfgenes are in
[Performance Model](performance.md#mgsa).

Differences from Ontologizer's MGSA: Ontologizer samples `alpha`, `beta`, and
`p` over a discrete grid alongside the term states, while this implementation
integrates them out against Beta priors. Expect similar rankings rather than
identical posteriors. On the `examples/ontologizer-ui` study and population
lists with current human GOA, all nine terms that the Ontologizer UI reports as
high-scoring rank in the top 12 here.

## Preranked GSEA

`genesets-rs gsea` implements the preranked statistic of Subramanian et al.,
[Gene set enrichment analysis: a knowledge-based approach for interpreting
genome-wide expression profiles](https://doi.org/10.1073/pnas.0506580102)
(2005). It is the right tool when the input is a quantity for every measured
gene (a fold change, a signed -log10 p-value, a t-statistic) rather than a list
chosen by a cutoff: no threshold is needed, and genes with larger scores count
for more.

Genes are sorted by score, highest first. For a target with `k` members among
`N` ranked genes, a running sum walks down the list, stepping up by
`|score|^weight / NR` at each member (`NR` sums that quantity over members) and
down by `1 / (N - k)` at every other gene. The enrichment score (ES) is the
running sum's largest deviation from zero: positive when members gather at the
top, negative at the bottom. `--weight 0` is the classic Kolmogorov-Smirnov
statistic, `--weight 1` (the default) is GSEA's weighted statistic. The leading
edge is the members up to the peak (positive ES) or from the trough on
(negative ES).

Significance uses gene-set permutation, as in GSEA preranked and fgsea's simple
mode. The null for size `k` is the ES of `--permutations` random `k`-gene
subsets of the same ranking, and every target with the same in-ranking size
shares that null, so cost scales with the number of distinct sizes. The
p-value is `(b + 1) / (m + 1)`, where `m` is the number of null draws with the
same sign as the observed ES and `b` the number at least as extreme. Its floor
is therefore about `2 / permutations`; raise `--permutations` when many terms
sit at the floor. NES divides ES by the mean of the same-sign null.

That floor is why `p_adjust` for KS is GSEA's FDR q-value rather than
Benjamini-Hochberg. At 1,000 permutations no term's p-value goes below about
0.002, so BH over the roughly 4,400 GO Biological Process terms in the 15-500
size window can never reach 0.05, however strong the signal. GSEA's FDR pools
the normalized null NES of every tested term and, for an observed NES, divides
the share of pooled null NES at least as extreme by the share of observed NES
at least as extreme, within each sign. Pooling resolves far below one term's
floor. The q-values are made monotone in |NES|, as BH is in p.

Each ranking and set size draws from its own seeded stream, so results do not
depend on `--threads`. Gene-set permutation does not preserve gene-gene
correlation, unlike GSEA's phenotype permutation, which needs the expression
matrix; this is the usual caveat for every preranked tool.

Checked against gseapy 1.3.1 `prerank` on a synthetic 3,000-gene ranking and 60
sets: ES agrees to within gseapy's output rounding (5e-7), NES correlates at
0.99998, and both flag the same 40 sets at nominal p < 0.05. On a random
20,000-gene ranking against 15,000 random sets, 4.9% of tested sets have
nominal p < 0.05 and none reach FDR q < 0.25. Rerun the gseapy check with
`uv run scripts/check_gsea_against_gseapy.py` after `cargo build --release`;
it also checks the FDR q-values against gseapy's.

### Wilcoxon rank-sum

`--statistic wilcoxon` tests members against non-members with the Wilcoxon
rank-sum (Mann-Whitney U) statistic: tie-corrected normal approximation with
continuity correction and a two-sided p-value, as scipy's
`mannwhitneyu(method="asymptotic")` and R's `wilcox.test(exact = FALSE)`
compute it. It needs no permutations, so p-values are not floored and BH
applies directly. `es` is the rank-biserial correlation, `2U / (k (N - k)) - 1`,
and `nes` the z-score. It weighs every member equally, where GSEA's weighted KS
concentrates on the extremes of the ranking; the leading edge reported for it
comes from the unweighted running sum.

### Using the ontology graph

By default (`--algorithm classic`) every term is tested on its own, so an
ancestor of a strongly enriched term is usually significant too, carried by
the child's genes. `--algorithm elim` and `--algorithm weight01` follow
Alexa, Rahnenführer & Lengauer, [Improved scoring of functional groups from
gene expression data by decorrelating GO graph
structure](https://doi.org/10.1093/bioinformatics/btl140) (2006), as the topGO
package implements them. Terms are scored deepest first, so every child is
final before its parent.

- `elim` scores a term without the genes of descendants whose raw p-value is
  at or below `--elim-cutoff` (default 0.01, topGO's default).
- `weight01` scores a term without the genes already removed below it, then,
  for each direct child at least as significant as the term, removes that
  child's genes and rescores, until no remaining child is. A parent that is
  more significant than all its children keeps their genes. "At least as
  significant" means a smaller p-value, or an equal one and an |NES| at least
  as large, which matters at the KS permutation floor.

In both, the removed genes leave the ranking for that term's test, as in
topGO's elim for rank statistics: the question is whether the term is enriched
among the genes not already explained below it. topGO's `weight` algorithm,
which down-weights rather than removes genes, is defined for Fisher's test
only, there and here.

Two caveats. These p-values are conditional on earlier results from the same
data, so they are not independent tests and their `p_adjust` is a heuristic
ranking aid rather than a calibrated FDR. And topGO 2.65's `.sigRatio.01`
(`R/topGOfunctions.R`) returns a ratio above 1 for every child, so topGO's own
`weight01` removes every child's genes regardless of significance, testing
each term on the genes annotated to it and to no child. The rule implemented
here is the conditional one topGO describes; expect it to differ from topGO's
`weight01` output.

For KS, a term with removed genes gets its own permutation null over its
reduced ranking, so graph-aware KS runs cost more than classic ones; Wilcoxon
stays analytic.

How these choices compare on the real GO graph is measured in
`evals/gsea_planted` (planted signal) and `evals/gsea_benchmark` (the curated
gold, with synthetic scores); their READMEs summarize the findings. In short:
both algorithms remove most redundant ancestors; elim loses a signal that
belongs to a general term (found in under 10% of planted runs against 66-85%
for classic), while weight01 with Wilcoxon keeps most of it; and both cut
recall of curator CORE terms roughly in half, because curators often place
core biology at levels these algorithms treat as redundant. `classic` stays
the default.

## Future Statistics

The code is structured so additional test families can be added without changing the data model:

- Benjamini-Hochberg FDR for `enrich` and `matrix` (`gsea` already reports it);
- multilevel (fgsea-style) KS p-values below the permutation floor;
- topology-aware parent-child methods;
- elim, weight and weight01 for Fisher `enrich`/`matrix` (`gsea` has elim and
  weight01).
