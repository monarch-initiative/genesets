# Ballouz, Pavlidis & Gillis 2017: multifunctionality case studies

Ballouz S, Pavlidis P, Gillis J. *Using predictive specificity to determine
when gene set analysis is biologically meaningful.* Nucleic Acids Res
2017;45(4):e20. PMID:28204549, PMCID:PMC5389513,
[doi:10.1093/nar/gkw957](https://doi.org/10.1093/nar/gkw957).

The paper shows that GO enrichment of several published gene lists is driven
by multifunctional (highly annotated) genes. The same case studies were reused
by the PAN-GO human functionome paper (Feuermann et al. 2025 Nature,
PMID:40011791, Supplementary Information) to argue that IBA-only annotations
avoid this bias. Here they are raw material for answer-keyed curation, which
records what the source authors judged meaningful **and** what they judged to
be artefacts (`category: false_association` / `nonspecific`).

## Files

`case_studies.gmt` holds one line per list:
`<name>\t<provenance>\t<genes...>`.

| set | genes | source of membership |
|---|---|---|
| GILMAN_2011_AUTISM_CNV | 72 | Gilman et al. 2011 *Neuron* (PMID:21658583) Table S2, copied from ai-gene-review `genesets/Gilman_autism_de_novo_CNV_genes_PMID21658583/genes.csv` |
| MANALO_2005_HYPOXIA_HIF1_UP | 202 | Manalo et al. 2005 *Blood* (PMID:15374877), ErmineJ quicklist from the paper's companion page |
| SCHMIDT_KASTNER_2012_SCHIZOPHRENIA | 77 | Schmidt-Kastner et al. 2012 *Mol Psychiatry* (PMID:22290124), ErmineJ quicklist |
| PARDO_2010_OCT4_INTERACTORS_MOUSE | 87 | Pardo et al. 2010 *Cell Stem Cell* (PMID:20362542), ErmineJ quicklist. **Mouse symbols**: these need ortholog mapping before human GOA enrichment |

The quicklists come from the companion page
<https://erminej.msl.ubc.ca/multifuncsupplement/>. They are not in the NAR or
PMC supplement, which is a methods PDF only.

**Gilman list.** The companion page's Gilman quicklist (`vitkup-cluster1`)
has only 44 genes. It lacks 5 of the 11 multifunctional genes the paper says
it removed (DLG1, CRHR1, DKK1, AXIN1, WNT3), so it is not the list the paper
analysed. The 72-gene Table S2 list used here contains all 11, and its size
matches the 72 genes reported by Feuermann et al.

**MSigDB lists.** The paper also tested 1,800 MSigDB v3.1 C2 lists (2013).
1,707 of them are in the current MSigDB C2 CGP sets fetched by the AIGR
change-set eval (`evals/aigr_changeset`). Their membership has drifted (median
Jaccard 0.74), so they are not copied here.

**Not obtained:** the "top 100 multifunctional genes" list the paper used to
compare methods. No file for it was found.

## What the authors judged meaningful vs artefactual

All quotes below were checked verbatim against the Europe PMC full text.

**Gilman autism CNV**
- **Artefact:** *learning or memory*. Quote: "heavy down-weighting of the
  'learning or memory' GO gene set". Four of its six genes are among the 11
  highly multifunctional genes (NRXN1, NLGN3, DLG4, CRHR1). Quote: "we prefer
  to think of it as non-robust and, most importantly, non-specific".
  Supplementary Figure 7 adds *single-organism behavior* and *positive
  regulation of signaling* as non-specific groups that lose significance.
- **Meaningful:** *neuron migration*, which is unaffected by the correction.
  Quote: it "might be even more relevant to ASD than learning and memory".

**Manalo hypoxia**
- **Meaningful:** 70 terms are significant before correction and 4 after:
  *peptidyl-proline modification*, *cellular response to hypoxia*,
  *collagen fibril organization* and *cellular response to oxygen levels*.
  Quote: these "strongly align with the themes the authors chose".
- **Artefact:** the other terms are implied artefacts, but the paper does not
  name them individually.

**Schmidt-Kastner schizophrenia**
- **Artefact:** GO enrichment results "related to stress and inflammation".
  Quote: a multifunctionality-corrected analysis "results in no gene sets
  meeting the significance criterion".
- **Meaningful:** none reaches significance. The paper's point is that "it is
  possible to construct a variety of narratives".

**Pardo Oct4 interactors**
- **Artefact:** named examples include *modification of symbiont morphology or
  physiology* and *ATP metabolic process*, plus terms relating to embryonic
  development.
- **Meaningful:** "chromatin remodeling and histone acetylation".

These judgements are the source authors' views, under their 2013–2014 GO and
their correction method. When curating, record them as evidence (with the
quote) on the relevant associations. Do not treat them as automatically
correct.
