# Model-organism gold standard: GO enrichment for mechanism

**Question:** When a gene list comes from a mouse, worm, fly, zebrafish, budding
yeast or fission yeast experiment, does GO enrichment against that species' own annotations recover
the curated biology, and in particular the *mechanism* behind the phenotype?

The rest of the gold standard is human. These 20 sets were chosen because their
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
| `GUO_TARGETS_OF_IRS1_AND_IRS2` | mouse | MSigDB M2:CGP | Hepatic insulin signalling sustains the sterol programme (cholesterol synthesis 14/35) |
| `ZHENG_IL22_SIGNALING_UP` | mouse | MSigDB M2:CGP | IL-22 drives epithelial antimicrobial defence; Reg3b/Reg3g missing (membership gap) |
| `WELCH_GATA1_TARGETS` | mouse | MSigDB M2:CGP | GATA-1 drives heme synthesis and erythroid maturation (largely confirmatory) |
| `CTD_ALDICARB_CELEGANS` | worm | CTD | Release machinery, plus Gq/Go, neuropeptide and dystrophin layers |
| `CTD_BT_CRYSTAL_PROTEIN_CELEGANS` | worm | CTD | p38/JNK toxin defence, GPI anchors, membrane repair (mostly GO gaps) |
| `CTD_PARAQUAT_CELEGANS` | worm | CTD | ETC, insulin/DAF-16, and sensory-cilium IFT |
| `CTD_PARAQUAT_DMEL` | fly | CTD | Detoxification, JNK signalling, dopamine/Parkinson's module |
| `CTD_TCDD_DRERIO` | zebrafish | CTD | AHR complex, fin regeneration, cardiac targets (a GO gap) |
| `YEAST_RLS_EXTENDED_MCCORMICK` | yeast | SGD / McCormick 2015 | Lifespan extension is 60S-specific (30 RPL vs 2 RPS genes), pointing to Gcn4 |
| `YEAST_RAPAMYCIN_SENSITIVE_XIE` | yeast | SGD / Xie 2005 | The vacuole as TORC1 hub: full EGO (Ragulator) complex, vacuolar sorting, autophagy |
| `YEAST_CADMIUM_SENSITIVE_RUOTOLO` | yeast | SGD / Ruotolo 2008 | Glutathione chelation, YCF1 vacuolar import, V-ATPase, ESCRT/retromer |
| `YEAST_NICKEL_RESISTANT_RUOTOLO` | yeast | SGD / Ruotolo 2008 | Losing ESCRT/MVB sorting protects against nickel |
| `YEAST_DTT_SENSITIVE_RAND` | yeast | SGD / Rand & Grant 2006 | Secretory pathway and UPR; V-ATPase as a nonspecific multidrug class |
| `POMBE_CADMIUM_SENSITIVE_KENNEDY` | fission yeast | PomBase / Kennedy 2008 | Phytochelatins, sulfate assimilation and ubiquinone, not budding yeast's GSH/YCF1 route |
| `POMBE_CLS_EXTENDED_ROMILA` | fission yeast | PomBase / Romila 2021 | PKA and TORC1 longevity genes present, but no GO term significant (hard case) |

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

The yeast sets are `LIT:GENETIC` screen convergences (one paper, one phenotype).
Each is the set of null mutants that SGD's literature curation records with that
phenotype from that paper (`phenotype_data.tab`, the "null" mutant type). Where SGD's
count differs from the paper's, the YAML says so (DTT: 174 curated versus 195
reported). Unlike the CTD sets, each has a single defining paper, so they carry
verbatim evidence snippets. The fission-yeast sets follow the same pattern using
PomBase's curated deletion phenotypes (FYPO terms) from one paper each, keeping
only screens whose PomBase count matches the paper's. That rule excluded two
candidates: bortezomib (255 curated vs 19 reported) and hydroxyurea (302 vs 52).

**Membership** for all 20 is snapshotted in
`curation/genesets/model_organism_members.gmt` (gene symbols as used by each
species' GAF; yeast uses SGD standard names, or the systematic ORF name when a
gene has no standard name). This is the analogue of `lit_members.gmt`, and makes the eval
reproducible if mygeneset.info changes.

## Pipeline

```bash
cargo build --release
uv run --project python/genesets-workflows --extra curation \
  python scripts/run_model_organism_eval.py --out-dir /tmp/mo_eval
```

For each taxon, the runner writes that taxon's sets to `queries.gmt` and runs
`scripts/prepare_go_eval.py` with the species GAF (`mgi`, `wb`, `fb`, `zfin`;
all evidence, NOT-filtered, symbol ids, a species-wide GOA background; `sgd` and
`pombase` for the yeasts). Non-human sets whose taxon has no GO Consortium GAF
(e.g. the Poeciliidae sets, which are mapped to human orthologs and scored in the
`goa_human` eval) are reported and skipped. It then
runs `genesets-rs run` (Bonferroni, `max_p_adjust 0.05`) and scores with
`scripts/score_method_vs_benchmark.py --taxon <NCBITaxon>`.

The runner passes `--gaf-db` with each resource's own DB column (`MGI`, `WB`, `FB`,
`ZFIN`, `SGD`, `PomBase`). `prepare_go_eval.py` keeps only `UniProtKB` rows by default, which
suits `goa_human` but would load almost nothing from a model-organism GAF: 311 of
MGI's ~707,000 rows, and none at all for SGD. Like the human default, this excludes
ComplexPortal complexes, and direct *protein binding* annotations are dropped.

The `--taxon` filter keeps human gold out of a species run and vice versa. Never
score a model-organism set against `goa_human`: its symbols would mostly fail to
match.

Bonferroni here corrects across every set × term in a run. Adjusted p-values
therefore depend on how many sets run together, so the curator notes quote
**nominal** p-values from this per-species run.

## Results (GO releases/2026-07-26; all-evidence GAFs)

| species | sets | recall_core | confirmatory | recall_confirm | mechanistic | recall_mechan |
|---|---|---|---|---|---|---|
| mouse (mgi) | 8 | 0.850 (17/20) | 19 | 0.684 | 13 | 0.692 |
| worm (wb) | 3 | 0.667 (6/9) | 8 | 0.625 | 10 | 0.200 |
| fly (fb) | 1 | 0.750 (3/4) | 4 | 0.750 | 2 | 0.500 |
| zebrafish (zfin) | 1 | 0.667 (2/3) | 2 | 1.000 | 3 | 0.333 |
| budding yeast (sgd) | 5 | 0.600 (6/10) | 7 | 0.571 | 13 | 0.538 |
| fission yeast (pombase) | 1 of 2 | 0.500 (2/4) | 2 | 0.500 | 3 | 0.333 |
| **pooled** | 19 of 20 | **0.720 (36/50)** | 42 | 0.667 | 44 | 0.477 |

`gap_recovered` is 0 in every species: no curator `annotation_gap` call was
contradicted by a tool run.

The scorer counts only sets with at least one significant hit. `POMBE_CLS_EXTENDED_ROMILA` has none, so it is not in the table and its three
unrecovered core and supporting terms do not lower recall. Read it as a failure
that recall does not show (finding 8).

### Findings

1. **When the carrier genes are annotated, enrichment recovers the mechanism.**
   Outside the yeasts, 28 of the 30 core terms predicted `annotation_supported`
   are recovered. The two misses are both small-term cases. One is the worm p38
   MAPK cascade, with 3 annotated carriers in an 8-gene term. The other is mouse
   gluconeogenesis in the IRS1/2 set (Pck1, G6pc1, Ppargc1a: 3/51). Each is a
   strong fraction with too few genes to survive genome-wide correction. Yeast shows the same failure mode at scale
   (finding 5). The mouse mechanistic terms (sterol synthesis behind POR loss and
   Egr2 loss, lipid homeostasis behind Atg16l1 loss, TOP-mRNA ribosome behind
   rapamycin) and the worm cilium/IFT and ETC terms behind paraquat resistance
   are all recovered.
2. **Outside the yeasts, six of the eight missed core terms are gaps:** five in
   species annotation, one in membership.
   These are biology the membership carries but species GOA does not: mouse
   acute-phase response (Saa1, Saa2, Hp present but unannotated, and again Hp,
   Trf, Lrg1 in the IL-22 set, so the gap is systematic in MGI); worm response
   to toxic substance (29 *hpo* Cry5B-hypersensitivity genes unannotated);
   response to paraquat in worm and fly (1 annotated member each, though most
   members entered CTD on paraquat phenotypes); and zebrafish heart
   development (tcf21, notch1b, aldh1a2 annotated elsewhere). Each names a
   concrete GO-annotation target. The paraquat cases show phenotype screens
   whose results never become process annotations. The membership gap is
   IL-22's Gram-negative defence: the paper's Reg3b and Reg3g are absent from the
   MSigDB set.
3. **Mechanistic recall is much worse in worm** (0.2) than in mouse (0.78) or
   yeast (0.54). Small
   signalling and modulatory modules are hard to recover: neuropeptide, Gq/Go
   and dystrophin layers in aldicarb resistance; p38/JNK, GPI anchors and
   membrane repair in Cry5B defence. They have few annotated genes, sit in large
   terms (GPCR signalling, 662 genes), or have their carriers annotated to
   sibling terms.
4. **Specific-sounding false positives come from moonlighting or screen-derived
   annotations.** Ribosomal proteins carry annotations to *negative regulation of
   myoblast fusion* and *translation at presynapse*. Both are significant in mouse
   embryonic fibroblasts and recorded as `false_association` negative controls.
   In yeast, ESCRT genes carry IMP annotations to *ATP export* from a single
   deletion screen of extracellular ATP (PMID:26585826). That makes "ATP export"
   the second-strongest hit for nickel resistance (15/17, p ~3e-29) and
   significant for cadmium and rapamycin too: screen phenotypes converted to
   process annotations leak into every ESCRT-rich list.
5. **In the yeasts, the specific mechanism often lives in a complete but tiny
   complex.** All six yeast misses are `annotation_supported`, with every carrier
   annotated. Fission yeast adds phytochelatin synthesis (pcs2) and sulfate
   assimilation (5/8) for cadmium. The four budding-yeast misses are the whole EGO/Ragulator complex for rapamycin (GTR1, GTR2,
   MEH1, SLM4: 4/6), TORC1 (TOR1, KOG1), glutathione synthesis for cadmium (GSH1,
   GSH2: 2/3), and the UPR for DTT (IRE1, HAC1). Genome-wide Bonferroni over
   thousands of terms cannot reward 100% coverage of a 2-6 gene complex. Instead
   enrichment reports the large surrounding compartment (vacuolar transport,
   vacuolar acidification, Golgi), which is partly a generic multidrug-sensitivity
   signal. Yeast is where small-term power matters most.
6. **The yeast lifespan phenotype has no GO term at all**: GO has obsoleted
   *replicative cell aging* and the other aging terms. The 60S-subunit hit is
   still recovered, so the mechanism is visible even though the phenotype is
   not.
7. **Two yeasts, two cadmium mechanisms, both visible to GO.** Budding yeast's
   cadmium screen recovers vacuolar acidification and ESCRT/retromer, with the
   specific glutathione/YCF1 route too small to score. Fission yeast's screen
   recovers ubiquinone biosynthesis (7/14, p ~1e-6) and response to metal ion.
   It also carries the phytochelatin route (pcs2, hmt1) and sulfate
   assimilation (5/8), again in terms too small to survive correction. The
   Kennedy paper itself notes how little the two screens overlap.
8. **A lifespan screen can carry its mechanism and still enrich for nothing.**
   The fission-yeast chronological-lifespan set contains the glucose/PKA (git3,
   git5, git11, pka1) and TORC1 (tco89, sck2) longevity genes, but not a single GO
   term passes Bonferroni. The pathways sit in 9- and 12-gene terms, and the
   remaining ~300 members are diffuse. It contrasts with the budding-yeast
   replicative-lifespan set, whose 60S signal is concentrated enough to recover.

### Caveats

- 20 sets (19 scored), 50 core terms: this is a pilot, not a benchmark-scale estimate.
- Because these sets were *selected* for mechanistic enrichment, 44 of 86
  scored insight-tagged terms are `mechanistic`, far above the ~10% corpus norm. Do
  not pool their mechanistic recall with the human corpus without saying so.
- CTD sets are aggregates; a few members may be artefacts of the underlying
  studies (noted per set).
- The guardrail is the same as `evals/iba_vs_benchmark`: these numbers measure
  the gold, and `recovery_status`/`insight` were set from GOA facts during
  curation, not refit to this run.
