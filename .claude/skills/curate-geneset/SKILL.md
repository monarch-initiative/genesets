---
name: curate-geneset
description: >-
  Use when curating or reviewing a GO-term interpretation of a non-GO gene set
  (MSigDB C8/C2/H, literature disease-activity, or GWAS/CRISPR genetics) for the
  curation/ gold standard in this repo. Covers grounding in the real membership,
  OLS term verification, the category / recovery_status / insight axes, evidence,
  series, the 4-gate validator, and folding the set into the eval.
---

# Curating a GO interpretation of a gene set

You are adding (or reviewing) one `curation/genesets/<SET>.yaml` — a curator's
GO-term interpretation of a non-GO gene set — to a hand-curated gold standard
that doubles as a precision/recall benchmark for `genesets-rs`. One YAML per set.

**Read first:** `curation/README.md` (the categories, recovery_status, insight,
sources, series, evidence — the authority) and the schema
`curation/schema/genesets_interpretation.yaml` (slots + enums). Match the house
style of an existing file, e.g. `curation/genesets/TRAVAGLINI_LUNG_CILIATED_CELL.yaml`.

## The one rule that governs everything: ground in the *actual* membership

**Before asserting any biology, read the set's real gene list.** Do not curate
from the cell type's textbook markers — curate from what is actually in the set.
This is what makes `recovery_status` honest and is the most common way agents go
wrong (e.g. a fetal "photoreceptor" set that lacks RHO/PDE6; a "chromaffin" set
that lacks TH/DBH; a "pro-B" set dominated by cell-cycle genes).

How to get the membership:
- **MSigDB sets** (`MSIGDB:<NAME>`, collections C8/C2/H/…): `_id` on
  mygeneset.info equals the MSigDB set name. Fetch with
  `scripts/fetch_mygeneset_query.py` (or query `_id:"<NAME>"` directly). The full
  list is also what gets folded into the eval.
- **Literature sets** (`LIT:` ids): the defining paper *is* the identity,
  membership, and evidence. For short, explicit lists (GWAS loci, CRISPR hits,
  causal-gene panels) capture the genes verbatim from the primary source,
  HGNC-normalize, and append to `curation/genesets/lit_members.gmt` (tab-sep:
  `<gene_set_name>\tcurated_literature_membership\tGENE1\tGENE2…`). The GMT
  `set_id` MUST equal the YAML `gene_set_name`.

## Per-term: the three orthogonal axes

For each association set `category`, `confidence`, `specificity`, `seed_source`,
`recovery_status`, and (on core/supporting terms) `insight`, plus a `curator_note`
and optional `evidence`. Aim for ~5–7 associations: the cell/state's defining
core terms (a process + a component), 1–2 supporting, at most one `nonspecific`
*only if* the membership clearly carries it (don't invent a housekeeping hit).

1. **`category` — the biology** (authoritative; NOT driven by annotation state).
   `core_process` / `core_component` / `supporting_process` /
   `marker_driven_plausible` / `nonspecific` / `false_association`. A term that is
   biologically core stays core even if no annotation supports it.

2. **`recovery_status` — grounded in the actual membership** (orthogonal to
   category). Name the specific present/absent genes in the `curator_note`:
   - `annotation_supported` — carrier genes are in the set; enrichment recovers it.
   - `annotation_gap` — relevant genes ARE in the set but GO annotation is too
     shallow. The gap is GO's (a GO-annotation curation target). `seed_source: curator_added`.
   - `membership_gap` — the carrier genes are NOT in the set (legacy/incomplete
     set, or a fetal signature lacking its mature effectors). The term still
     belongs in a complete description. `seed_source: curator_added`.

3. **`insight` — interpretive value** (only on core/supporting). The tightened rule:
   - `confirmatory` — entailed by the set's construction/identity (the default for
     a marker signature; ALL generic hallmarks — proliferation, apoptosis, generic
     PI3K/MAPK, known disease OXPHOS, expected immune terms).
   - `mechanistic` — a SPECIFIC, non-obvious process not entailed by construction;
     a genuine enrichment insight. Use sparingly (~10% of terms corpus-wide). When
     in doubt, it's confirmatory.

## Anti-hallucination discipline (non-negotiable)

- **OLS-verify every ontology id+label** (GO/CL/UBERON/MONDO/…) before writing it.
  Use the OLS MCP (`mcp__claude_ai_OLS__search` / `fetch`). Confirm the id resolves
  to the *exact* label and is **not obsolete**. Never invent ids; if unsure, search
  by label and use the returned id.
- **Reject obsolete terms.** The validator's 4th gate sweeps OAK `obsoletes()` and
  fails the build. The `sqlite:obo` builds retain labels for obsolete classes, so
  an obsolete id paired with its old label passes id+label validation but fails
  here (e.g. `GO:0062023` collagen-containing ECM → `GO:0031012`; `GO:0050663`
  cytokine secretion → `GO:0032635`).
- **Evidence snippets are verbatim or absent.** Include an `evidence` item with a
  `snippet` ONLY if you retrieved the exact text from the cited paper (PubMed MCP)
  and copied it character-for-character — the reference-validator substring-checks
  it and fails on any mismatch. Otherwise use `curator_note` only. Never fabricate.
  Any item with a `snippet` must also carry a `reference`.
- **Cite the right paper.** Verify PMIDs. (The Descartes fetal atlas is
  `PMID:33184181`, not `32848094` — a recurring trap.)

## Pairs / series

Link contrasting poles of one axis with `series` (a shared `SERIES:<NAME>` id) and
`series_role` (this set's pole, free text). Add the field to *both/all* poles
(edit the partner files too). The eval checks that opposite poles resolve to
contrasting GO interpretations — e.g. `SERIES:PANCREATIC_ISLET` alpha/beta/delta,
`SERIES:MYELINATING_GLIA` central/peripheral.

## Validate, then regenerate the manifest

```bash
# 4 gates: linkml structural, term id+label, reference snippet, obsolescence
uv run --project python/genesets-workflows --extra curation \
  genesets-workflows curate validate curation/genesets/<SET>.yaml   # exit 0 = pass
```
Iterate until exit 0. Validate the whole corpus before a PR with `just curate-validate`.

Regenerate `curation/genesets/manifest.tsv` after a batch (columns:
`gene_set_id, collection, context_term, context_type, context_label, series,
series_role`) by looping every YAML through `model.load_interpretation`. Run
`just curate-test` (29 unit + 5 doctest).

## Fold the set into the eval (so it becomes scoreable)

A set is only evaluable if its membership is in the eval's `queries.gmt`. MSigDB
sets come from the mygeneset base; `LIT:` sets come from `lit_members.gmt`. After
adding sets, rebuild `queries.gmt` and re-run the **evaluate-enrichment** skill.

## The guardrail (read `evals/iba_vs_benchmark/README.md`)

The eval *measures* the gold; it must never *refit* it. `category` is the scored
truth. `recovery_status`/`insight` are curator judgments adjudicated against GOA
facts during curation — never auto-updated to match whatever a tool happened to
recover (that would make recall-vs-gold circular). When the eval surfaces a
disagreement (a `gap_recovered` term), it is a review item, decided on the merits,
not an automatic relabel.
