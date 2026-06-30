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

## Gene-Set-Informed Factor Models

A separate paradigm consumes gene sets without doing over-representation analysis
at all. Instead of asking "is this set enriched in my ranked or thresholded
gene list?", these methods treat curated sets as *priors* on the latent
structure of an expression matrix and let the model decide, per sample or per
cell, which annotated programs are active.

f-scLVM (Buettner et al., 2017), implemented in the **slalom** package for R and
Python, is the reference example. It is a sparse Bayesian factor analysis model
for single-cell RNA-seq in which each annotated factor is initialized from a gene
set (typically MSigDB or Reactome pathways): the set membership seeds the factor's
loading mask. The model then, jointly across all factors:

- estimates a per-factor relevance so it can switch annotated factors on or off
  rather than forcing every pathway to explain variance;
- refines the gene-set annotation itself, adding genes that load on a factor but
  were not in the input set and down-weighting members that do not;
- infers additional *unannotated* sparse and dense factors to absorb residual
  structure (technical and biological) that no input set explains.

The output is interpretable per-cell factor activations attached to named
pathways, which is why it is used to identify and annotate cell subpopulations
rather than to rank pathways for a single contrast.

expiMap (Lotfollahi et al., 2023), implemented in the **scArches** ecosystem, is
the deep-learning successor to the same idea. It is a variational autoencoder
whose latent nodes are gene programs: a binary mask ties the decoder so that each
latent dimension reconstructs only the genes of its annotated set, making the
otherwise opaque latent space directly interpretable. It carries over every move
f-scLVM makes — it learns a per-cell activity for each program, soft-prunes
uninformative programs (group-lasso regularization on the decoder, the relevance
analog), and refines membership by allowing the mask to admit a few extra genes
per program. Its headline addition is a set of **add-on de novo learnable
nodes**: unconstrained latent dimensions that capture biology the annotations
miss, the neural-network counterpart of f-scLVM's unannotated factors. Because it
is built on scArches, its distinctive use case is *reference mapping* — projecting
a new query dataset onto an existing atlas's gene-program latent space and reading
off which programs differ, including in perturbation and disease-response
settings.

The trade-off between the two mirrors linear-vs-nonlinear generally: f-scLVM is a
linear factor model that is fast, fully Bayesian, and easy to reason about;
expiMap is a nonlinear VAE that scales to atlas-sized data and supports transfer
learning, at the cost of interpretability being mediated by the architecture
rather than read directly off loadings.

This matters to a gene-set engine for two reasons. First, it is a downstream
*consumer* of exactly the libraries this project scores and curates, so set
quality propagates into factor interpretability. Second, the "refine the
annotation" step is the same concern the curated corpus encodes by hand: the
nominal membership of a named set and the genes that actually behave like the set
are not identical. Both f-scLVM and expiMap learn that gap from data (membership
refinement); the curation gold standard records it from expert judgment
(`recovery_status` / `membership_gap`). The two are complementary views of
annotation noise, not competitors.

GSEA and ssGSEA (single-sample GSEA) sit between ORA and factor models: they are
rank-based set-scoring methods that, in the ssGSEA case, also produce a per-sample
per-set activity score. They share the "score named programs per sample" goal with
f-scLVM but do not jointly model factors, switch sets off, or refine membership.

### Connection To The Curated Corpus

The curated gold standard in `curation/` is, structurally, a hand-built version
of what these models estimate. A single non-GO gene set is not curated to one
label; it is decomposed into *multiple* GO programs, each with a curator judgment.
`HALLMARK_INTERFERON_GAMMA_RESPONSE`, for example, resolves into a response-to-
type-II-interferon program, an antigen-processing/MHC-class-I program, an
antiviral-defense program, and a `nonspecific` translation residual. That is a
factor decomposition of the set, recorded with evidence rather than fit from an
expression matrix, and the schema axes line up with the model concepts:

- multiple `associations` per set ≈ the several latent program nodes one set seeds;
- `category` (`core_process` … `nonspecific` / `false_association`) ≈ per-factor
  relevance / group-lasso pruning — the `nonspecific` translation term is the
  factor f-scLVM would drive toward zero relevance;
- `recovery_status: membership_gap` ≈ data-driven **membership refinement**:
  expiMap's reported addition of nine B-cell markers to its predefined BCR program
  is exactly a `membership_gap` correction, learned from data instead of curated;
- `series` / contrasting poles ≈ a program's activation contrasting across cell
  states.

Two consequences. First, the corpus is natural **supervision and a benchmark for
the membership-refinement step specifically**: the models propose "gene X belongs
in program P"; the gold standard holds independent, cited judgments of that same
claim. Second, the alignment is only with `membership_gap` — the models refine set
membership but never touch the ontology, so `annotation_gap` (GO too shallow) stays
orthogonal, just as it does for the enrichment tools above.

The overlap is concrete, not aspirational: f-scLVM seeds factors from MSigDB
Hallmark and Reactome and headlines a G2/M-checkpoint cell-cycle factor; expiMap
trains on Reactome plus PanglaoDB marker sets and headlines interferon programs.
The corpus already curates that family — `HALLMARK_G2M_CHECKPOINT`,
`HALLMARK_INTERFERON_GAMMA_RESPONSE`, and the `HAY_BONE_MARROW_*` / `DESCARTES_*`
cell-type marker sets — so the priors these models consume are the same objects
this project scores. Note the scope line, though: this repo curates gene sets and
their GO interpretations, not expression matrices. The in-scope way to bring a
*signature* from one of these papers into the corpus is the existing
`LIT:DISEASE_ACTIVITY` pattern (a derived gene list with the paper as identity,
membership, and evidence), not hosting the underlying single-cell data.

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
- f-scLVM / slalom paper: <https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1334-8>
- expiMap paper: <https://www.nature.com/articles/s41556-022-01072-x>
- Enrichr paper: <https://bmcbioinformatics.biomedcentral.com/articles/10.1186/1471-2105-14-128>
- GO Consortium enrichment page: <https://geneontology.org/docs/go-enrichment-analysis/>
- MyGeneset.info docs: <https://docs.mygeneset.info/>
- OAK CLI docs: <https://incatools.github.io/ontology-access-kit/cli.html>
- Horned-OWL docs: <https://docs.rs/horned-owl/latest/horned_owl/>
