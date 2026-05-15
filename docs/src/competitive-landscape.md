# Competitive Landscape

The goal is not to be the broadest enrichment product. It is to be a very fast, reproducible core engine with an eval framework.

## Existing Tools And What They Optimize

GOATOOLS is a Python library and CLI for GO analysis. It parses OBO, GAF, GPAD, NCBI gene2go, and related formats; it propagates counts to parent terms by default; and it uses Fisher exact tests with many correction options.

g:Profiler and `gprofiler2` provide a maintained web/API/R workflow across GO, pathways, disease annotations, and identifier namespaces. They use hypergeometric enrichment and offer g:SCS, Bonferroni, and FDR correction choices.

clusterProfiler is the dominant R/Bioconductor ecosystem tool for ORA and ranked-list workflows. It has GO-specific helpers and a universal `enricher` interface using user-provided `TERM2GENE` mappings.

topGO is important because it exposes topology-aware GO methods such as `classic`, `elim`, `weight`, `weight01`, `lea`, and `parentchild` with multiple statistics.

Ontologizer is important for parent-child and model-set ideas. Its papers frame the issue that GO terms are dependent through inheritance, so standard term-for-term Fisher tests can be biologically redundant.

Enrichr is a popular web/API ecosystem with many gene set libraries, fast interactive workflows, and p-value, z-score, combined-score, and adjusted p-value outputs.

PANTHER powers the GO Consortium enrichment service and is a strong reference point for user expectations around GO term enrichment and background selection.

MyGeneset.info is not a direct competitor. It is a useful public gene-set source and API for eval inputs because it exposes GO, MSigDB, Disease Ontology/OMIM, WikiPathways, Reactome, KEGG, and species filters.

## Our Wedge

The core differentiators should be:

- ontology-neutral input tables rather than GO-specific parsing in the hot path;
- fast dense bitset scoring for repeated jobs;
- first-class N x N term/set matrix mode;
- explicit background semantics;
- reproducible prep metadata;
- benchmark and correctness evals as part of the repo.

## What Not To Compete On Yet

We should not try to match all visualization, identifier conversion, web UI, or database coverage of mature ecosystems in the core crate. Those are wrapper-layer concerns.

The MVP should instead make it easy to prove:

- this input produced this table;
- this table used this closure policy;
- these p-values match the reference implementation;
- this workload scored faster than comparable tools.

## Research Sources

- GOATOOLS project page: <https://pypi.org/project/goatools/>
- gprofiler2 vignette: <https://rdrr.io/cran/gprofiler2/f/inst/doc/gprofiler2.Rmd>
- topGO manual: <https://bioconductor.org/packages/devel/bioc/vignettes/topGO/inst/doc/topGO_manual.html>
- Ontologizer paper: <https://academic.oup.com/bioinformatics/article/24/14/1650/182451>
- Enrichr paper: <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128>
- GO Consortium enrichment page: <https://geneontology.org/docs/go-enrichment-analysis/>
- MyGeneset.info docs: <https://docs.mygeneset.info/>
- OAK CLI docs: <https://incatools.github.io/ontology-access-kit/cli.html>
- Horned-OWL docs: <https://docs.rs/horned-owl/latest/horned_owl/>
