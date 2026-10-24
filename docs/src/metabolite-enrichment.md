# Metabolite Enrichment

GO annotates gene products, not chemicals. A metabolomics result reaches GO
only through an intermediary, and the choice of intermediary decides what the
enrichment can say. `prepare-chem2go` materializes three such bridges as
ordinary genesets-rs annotation tables whose "genes" are canonical ChEBI ids,
so the engine runs unchanged.

| Bridge | Path | Reaches | Species |
|---|---|---|---|
| `logical_definition` | GO-plus axioms naming a ChEBI class (has_primary_input/output, has_participant, regulates_levels_of) | Chemical-named processes and transporters: *fatty acid metabolic process*, *citrate transport* | Independent |
| `rhea` | Rhea reaction participants → `rhea2go` | Catalytic molecular functions | Independent |
| `enzyme` (optional) | Rhea → UniProt enzyme → that enzyme's GAF annotations (BP by default) | Pathway processes: *tricarboxylic acid cycle*, *urea cycle* | One organism's GAF |

The bridges are kept separate so they can be compared the way GOA evidence
variants are (`annotations.<bridge>.tsv`, plus `all`, plus
`species_independent` when the enzyme bridge is on).

## Prepare

```bash
just chem2go-prep           # human GAF for the enzyme bridge
# or, species-independent bridges only:
uv run --project python/genesets-workflows genesets-workflows prepare-chem2go \
  --out-dir evals/chem2go/generated/current
```

Downloads (about 200 MB, cached under `downloads/`): GO-plus JSON, ChEBI lite,
ChEBI cross-references, the Rhea flat file, `rhea2go`, Rhea's pH 7.3 mapping,
and with `--gaf`/`--gaf-url` also `rhea2uniprot_sprot` and the GAF.

## Normalize, then enrich

```bash
genesets-workflows normalize-chemicals changed.txt  --prep-dir $PREP --out query.txt
genesets-workflows normalize-chemicals measured.txt --prep-dir $PREP --out background.txt
genesets-rs enrich --annotations $PREP/annotations.all.tsv \
  --terms $PREP/terms.tsv --closure $PREP/closure.tsv \
  --sample query.txt --background background.txt --min-overlap 2
```

`normalize-chemicals` accepts ChEBI in any spelling, HMDB (5- or 7-digit),
KEGG COMPOUND, LIPID MAPS, MetaCyc, DrugBank and CAS. Its report says how
each id mapped (`canonical`, `merged`, `secondary`, `xref`,
`xref_ambiguous`, ...) and whether the result carries any annotation.
`--extra-xrefs` takes a two-column TSV for ids ChEBI does not cross-reference.

## Design decisions

**Collapse what metabolomics cannot distinguish.** Conjugate acids/bases,
tautomers and Rhea's pH 7.3 mapping are merged, and the group is represented by
its pH 7.3 major microspecies, the form GO and Rhea use. Citric acid, citrate(3-)
and `HMDB0000094` all become `CHEBI:16947`. Enantiomers are *not* merged.

**Roll chemicals up the ChEBI hierarchy, but not through near-root classes.**
A chemical reaches a GO term whose axiom names it *or any ChEBI `is_a`
ancestor*, so palmitate reaches fatty acid metabolism. Roll-up stops at classes
with more than `--max-rollup-class-size` (default 2500) descendants; without
the limit "chemical entity" links every chemical to dozens of transporter and
homeostasis terms. Broad GO terms are still reached through the GO closure.

**Drop currency metabolites as participants.** Water, H+, ATP, NAD(P)H, CoA,
SAM and common ions are excluded wherever a chemical is one participant among
many (Rhea, `has_participant`). Primary inputs/outputs keep them, so
*ATP metabolic process* still reaches ATP. `currency.tsv` lists them with their
Rhea degree; `--keep-currency` disables this.

**Default universe = measurable metabolites.** Without `--universe`, only
chemicals that a bridge links directly, Rhea participants, and chemicals with
an HMDB, KEGG, LIPID MAPS or MetaCyc cross-reference are annotated (about 49k
of ChEBI's 220k classes; tripeptides and the like are left out).

**Link kinds.** `logical_definition` keeps metabolism, transport, catalysis and
homeostasis terms by default. *Binding*, *response to* and *signaling* terms
are dropped (`--link-kinds`): "response to glucose" describes the genes that
respond, not what happens to glucose.

## The background matters

The background must be the metabolites your platform measured, normalized with
the same command. The annotated chemicals are not a valid background:
untargeted platforms detect a biased fraction of metabolism, and using all
annotated chemicals inflates significance for well-measured pathways (see
Wieder et al. 2021, *PLoS Comput Biol* 17:e1009105).

## What it recovers

On the smoke queries in `evals/chem2go/` (TCA cycle and urea cycle
intermediates, as mixed HMDB/KEGG/ChEBI ids):

- `logical_definition` and `rhea` return transporters and catalytic activities
  (C4-dicarboxylate transport, arginase activity), not the cycles themselves.
  GO rarely links pathway processes to chemicals: *tricarboxylic acid cycle*
  lists only citrate as a participant.
- `enzyme` puts *tricarboxylic acid cycle* and *urea cycle* in its top two,
  with some noise from multifunctional enzymes' other BP annotations
  (e.g. DNA repair terms on TCA enzymes).

So the species-independent bridges answer "which chemistry is this", and the
enzyme bridge answers "which pathway". Use `all` when you want both.

## Known gaps

- ChEBI's cross-references are uneven: `HMDB0000052` (argininosuccinate) is
  attached only to an unconnected 2-star stub, and carbamoyl phosphate's HMDB
  id is missing. `--extra-xrefs` patches these; a UniChem or MetaNetX export is
  the general fix.
- Lipid species reported at sum composition ("PC 34:1") need mapping to a
  ChEBI class first; once mapped, roll-up handles them.
- The enzyme bridge inherits GAF annotation biases, and one metabolite fans out
  to every pathway of every enzyme that touches it. Counting stays per
  metabolite, but expect broad terms.
- There is no curated metabolite benchmark yet. Natural answer keys: inborn
  errors of metabolism (known pathway, known accumulating metabolites),
  leave-one-out over Reactome/Rhea pathways, and published MetaboLights studies.
