# Statistics

The MVP performs standard one-sided over-representation analysis.

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

## Future Statistics

The code is structured so additional test families can be added without changing the data model:

- Benjamini-Hochberg FDR;
- ranked-list methods;
- topology-aware parent-child methods;
- elim and weight-style algorithms;
- model-set approaches inspired by Ontologizer.
