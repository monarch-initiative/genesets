# Post-Processing

GO enrichment routinely returns many related terms. The scoring engine should stay simple and complete; reduction should be a composable post-processing layer that consumes enrichment TSV plus closure/annotation metadata.

## Standard Families

Significance cutoff is the first reduction step. PANTHER-style APIs expose a cutoff parameter; when FDR correction is used, the cutoff is applied to FDR, otherwise to p-value. Our eval configs mirror this with `max_p_adjust: 0.05`.

GO slims are a curated reduction strategy. Instead of pruning enriched terms after scoring, project annotations or results onto a smaller subset of broad terms. This is useful for overviews, less useful when the user wants precise mechanistic terms.

Slim-descendant filters are a related report-ranking strategy. Instead of
projecting results onto slim terms, keep candidate terms that are themselves in
a slim or descend from a slim term, often with a target-size guard. This is
useful for "largest changed term" summaries where root-like terms such as broad
molecular function categories would otherwise dominate.

Antislim filters invert that idea: exclude a slim region and its descendants
from a report. This is useful when a known broad ontology branch is technically
correct but distracts from the biological question being reviewed.

Ancestor-descendant pruning uses the closure graph. A simple policy is: for each query, if a term has a significant descendant with equal or better adjusted p-value, hide the ancestor. This is easy to explain, but it can over-prune when a broad parent captures a coherent process and many children are sparse.

Gene-overlap pruning compares term gene sets. If two enriched terms share nearly the same query-overlap genes, keep one representative. Representative selection can rank by adjusted p-value, specificity, overlap size, or information content.

GO Trimming-style approaches remove redundant parent terms based on graph paths and result statistics after enrichment has already been run.

REVIGO-style approaches cluster terms by semantic similarity and keep representatives. This is more flexible than strict ancestor-descendant pruning because related terms need not sit on the same direct path.

Topology-aware scoring methods such as topGO `elim`, `weight`, `weight01`, and `parentchild`, or Ontologizer-style parent-child/model-set methods, are not just post-processing. They change the p-values by accounting for the GO graph during scoring. We should support them later as alternate statistical models, not as TSV filters.

## Proposed Composable Reducers

Reducers should take a full result table and emit:

- a reduced table;
- a mapping from hidden term to representative term;
- a reason code;
- parameters and input file digests in metadata.

Candidate reducer modes:

| Mode | Keep | Hide |
| --- | --- | --- |
| `ancestor-dominated` | more specific descendant | ancestor when descendant has p-adjust <= ancestor p-adjust times tolerance |
| `ancestor-dominated-strict` | more significant term | ancestor only when descendant has equal or better p-adjust |
| `overlap-jaccard` | best ranked representative | terms whose overlap gene Jaccard exceeds threshold |
| `semantic-cluster` | cluster representative | terms in same semantic-similarity cluster |
| `go-slim-project` | slim term | non-slim terms after projection |
| `slim-descendant-include` | terms in or below selected slim terms | terms outside the selected slim scope |
| `slim-descendant-exclude` | terms outside selected slim terms and descendants | terms in or below selected slim terms |

The tolerance is important. Sometimes a general parent has a slightly better p-value because it adds many relevant genes. A practical rule is not binary dominance, but dominance with a margin:

```text
hide ancestor A if descendant D is significant and
  p_adjust(D) <= p_adjust(A) * tolerance
```

With `tolerance = 1.0`, the descendant must be at least as significant. With `tolerance = 2.0`, a slightly worse but much more specific descendant can still represent the signal.

## PANTHER Behavior

PANTHER's public service exposes Fisher or binomial tests, FDR/Bonferroni/none correction, and a cutoff parameter. Public documentation emphasizes cutoff-filtered result tables and GO-slim datasets. I did not find evidence that PANTHER applies a dynamic ancestor-descendant redundancy-pruning pass to complete GO over-representation outputs. For our purposes, treat PANTHER as a reference for cutoff filtering and GO-slim-style summaries, not as a specific redundancy-pruning algorithm.

## Recommended Default

For eval tables, keep the unpruned significant TSV as the canonical result:

```yaml
max_p_adjust: 0.05
post_processing: none
```

Then create optional reduced views:

```yaml
post_processing:
  method: ancestor-dominated
  p_adjust_tolerance: 2.0
  prefer: specificity
```

This avoids losing terms before diffing or debugging.

## Sources

- PANTHER API parameters via rbioapi docs: <https://rbioapi.moosa-r.com/reference/rba_panther_enrich.html>
- GO subset guide: <https://geneontology.org/docs/go-subset-guide/>
- GO Trimming paper: <https://pmc.ncbi.nlm.nih.gov/articles/PMC3160396/>
- REVIGO paper: <https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0021800>
- topGO manual: <https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO_manual.html>
