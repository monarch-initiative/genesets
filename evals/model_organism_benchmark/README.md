# Model-organism gold standard: GO enrichment for mechanism

**Question:** When a gene list comes from a mouse, worm, fly or zebrafish
experiment, does GO enrichment against that species' own annotations recover
the curated biology, and in particular the *mechanism* behind the phenotype?

The rest of the gold standard is human. These 10 sets were chosen because their
enrichment should point to a mechanism rather than restate how the set was
built. For example, deleting a generic P450 reductase turns out to act through
cholesterol synthesis; an autophagy mutant reveals a lipid programme; and
aldicarb resistance resolves into exocytosis, Gq/Go and neuropeptide layers.
Each set is grounded in its real membership.

## The sets

| set | species | source | what enrichment reveals |
|---|---|---|---|
| `SCHMIDT_POR_TARGETS_IN_LIMB_BUD_UP` | mouse | MSigDB M2:CGP | POR loss acts through sterol synthesis and SREBP-2 feedback |
| `LE_EGR2_TARGETS_DN` | mouse | MSigDB M2:CGP | Egr2 drives myelin via cholesterol and galactolipid synthesis |
| `CADWELL_ATG16L1_TARGETS_UP` | mouse | MSigDB M2:CGP | Autophagy loss in Paneth cells gives a lipid and lipoprotein gain of function |
| `BILANGES_RAPAMYCIN_SENSITIVE_GENES` | mouse | MSigDB M2:CGP | mTORC1 acts through TOP-mRNA (ribosomal protein) translation |
| `LIAN_LIPA_TARGETS_6M` | mouse | MSigDB M2:CGP | Lysosomal storage, collagen-degrading MMPs, emphysema |
| `CTD_ALDICARB_CELEGANS` | worm | CTD | Release machinery, plus Gq/Go, neuropeptide and dystrophin layers |
| `CTD_BT_CRYSTAL_PROTEIN_CELEGANS` | worm | CTD | p38/JNK toxin defence, GPI anchors, membrane repair (mostly GO gaps) |
| `CTD_PARAQUAT_CELEGANS` | worm | CTD | ETC, insulin/DAF-16, and sensory-cilium IFT |
| `CTD_PARAQUAT_DMEL` | fly | CTD | Detoxification, JNK signalling, dopamine/Parkinson's module |
| `CTD_TCDD_DRERIO` | zebrafish | CTD | AHR complex, fin regeneration, cardiac targets (a GO gap) |

**Sources.** The mouse sets are native-mouse MSigDB M2 chemical/genetic
perturbation sets (`MSIGDB:` ids, `msigdb_release: v2025.1.Mm`). The worm, fly
and zebrafish sets are Comparative Toxicogenomics Database chemical-gene
interaction sets as served by mygeneset.info (`MYGENESET:<MeSH>_<taxid>`,
`collection: CTD:chem_gene_ixns`). CTD sets pool heterogeneous literature
(mutant sensitivity phenotypes, expression changes, proteomics) of mixed
direction, so they are `direction: NA` and carry no evidence snippets: there is
no single defining paper. Their curator notes say which members come from which
kind of study where that matters, for example the abundant-protein proteomic
signal in the TCDD set, recorded as `nonspecific`.

**Membership** for all 10 is snapshotted in
`curation/genesets/model_organism_members.gmt` (gene symbols as used by each
species' GAF). This is the analogue of `lit_members.gmt`, and makes the eval
reproducible if mygeneset.info changes.

## Pipeline

```bash
cargo build --release
uv run --project python/genesets-workflows --extra curation \
  python scripts/run_model_organism_eval.py --out-dir /tmp/mo_eval
```

For each taxon, the runner writes that taxon's sets to `queries.gmt` and runs
`scripts/prepare_go_eval.py` with the species GAF (`mgi`, `wb`, `fb`, `zfin`;
all evidence, NOT-filtered, symbol ids, a species-wide GOA background). It then
runs `genesets-rs run` (Bonferroni, `max_p_adjust 0.05`) and scores with
`scripts/score_method_vs_benchmark.py --taxon <NCBITaxon>`.

The `--taxon` filter keeps human gold out of a species run and vice versa. Never
score a model-organism set against `goa_human`: its symbols would mostly fail to
match.

Bonferroni here corrects across every set × term in a run. Adjusted p-values
therefore depend on how many sets run together, so the curator notes quote
**nominal** p-values from this per-species run.

## Results (GO current, 2026-11; all-evidence GAFs)

| species | sets | recall_core | confirmatory | recall_confirm | mechanistic | recall_mechan |
|---|---|---|---|---|---|---|
| mouse (mgi) | 5 | 0.923 (12/13) | 11 | 0.636 | 9 | 0.778 |
| worm (wb) | 3 | 0.667 (6/9) | 8 | 0.625 | 10 | 0.200 |
| fly (fb) | 1 | 0.750 (3/4) | 4 | 0.750 | 2 | 0.500 |
| zebrafish (zfin) | 1 | 0.667 (2/3) | 2 | 1.000 | 3 | 0.333 |
| **pooled** | 10 | **0.793 (23/29)** | 25 | 0.680 | 24 | 0.458 |

`gap_recovered` is 0 in every species: no curator `annotation_gap` call was
contradicted by a tool run.

### Findings

1. **When the carrier genes are annotated, enrichment recovers the mechanism.**
   Of the 24 core terms predicted `annotation_supported`, 23 are recovered. The
   one miss is the worm p38 MAPK cascade, which has 3 annotated carriers in an
   8-gene term: a strong fraction, but too few genes to survive genome-wide
   correction. The mouse mechanistic terms (sterol synthesis behind POR loss and
   Egr2 loss, lipid homeostasis behind Atg16l1 loss, TOP-mRNA ribosome behind
   rapamycin) and the worm cilium/IFT and ETC terms behind paraquat resistance
   are all recovered.
2. **Five of the six missed core terms are model-organism annotation gaps.**
   These are biology the membership carries but species GOA does not: mouse
   acute-phase response (Saa1, Saa2, Hp present but unannotated); worm response
   to toxic substance (29 *hpo* Cry5B-hypersensitivity genes unannotated);
   response to paraquat in worm and fly (1 annotated member each, though most
   members entered CTD on paraquat phenotypes); and zebrafish heart
   development (tcf21, notch1b, aldh1a2 annotated elsewhere). Each names a
   concrete GO-annotation target. The paraquat cases show phenotype screens
   whose results never become process annotations.
3. **Mechanistic recall is much worse outside mouse** (worm 0.2). Small
   signalling and modulatory modules are hard to recover: neuropeptide, Gq/Go
   and dystrophin layers in aldicarb resistance; p38/JNK, GPI anchors and
   membrane repair in Cry5B defence. They have few annotated genes, sit in large
   terms (GPCR signalling, 662 genes), or have their carriers annotated to
   sibling terms.
4. **Specific-sounding false positives come from moonlighting annotations.**
   Ribosomal proteins carry annotations to *negative regulation of myoblast
   fusion* and *translation at presynapse*. Both are significant in mouse
   embryonic fibroblasts and recorded as `false_association` negative controls.

### Caveats

- 10 sets, 29 core terms: this is a pilot, not a benchmark-scale estimate.
- Because these sets were *selected* for mechanistic enrichment, 24 of 49
  insight-tagged terms are `mechanistic`, far above the ~10% corpus norm. Do
  not pool their mechanistic recall with the human corpus without saying so.
- CTD sets are aggregates; a few members may be artefacts of the underlying
  studies (noted per set).
- The guardrail is the same as `evals/iba_vs_benchmark`: these numbers measure
  the gold, and `recovery_status`/`insight` were set from GOA facts during
  curation, not refit to this run.
