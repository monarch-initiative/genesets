# Semi-automated cGSA curation

Curating a cGSA-benchmark set has a **deterministic** half (the gene list) and a
**judgment** half (GO-term roles and context). `draft_cgsa.py` does the first and
hands you the second.

## What the tool does (deterministic)

For each benchmark set id it:

1. pulls the DEG membership + species + contexts from `evals/cgsa/data/`;
2. runs the repo's real enrichment kernel (`genesets-rs`) against **species-matched
   GOA** (human `goa_human`, mouse `MGI`), so every candidate GO term arrives with
   its true **overlap (carrier) genes** and Bonferroni stats — no hand-asserted
   carriers, no guessing which genes are annotated to what;
3. writes a **draft** interpretation — one association per significant term, with
   `seed_source: enrichment_recovered`, `recovery_status: annotation_supported`,
   the carriers in `enrichment_stats.overlap_genes`, and **`category` left unset**;
4. **screens for GT/membership mismatch**: it grounds the benchmark's own
   ground-truth pathway names to GO (via the grounding cache) and reports how many
   are actually recovered by enrichment. A low fraction flags the set (the
   `PANC_IMMUNE` failure mode — an epithelial DEG list with an all-immune ground
   truth) so you don't curate biology the membership doesn't support;
5. auto-suggests a disease context (top MONDO hit for the keywords) tagged
   `AUTO-SUGGESTED … verify/replace` — a candidate, never silently accepted.

```bash
python evals/cgsa/draft/draft_cgsa.py \
  --ids PMC9750880-1,PMC8983726-1 \
  --out-dir curation/cgsa/drafts --cache-dir /tmp/cgsa_annot
# -> curation/cgsa/drafts/CGSA_<id>.draft.yaml  +  screening_report.tsv
```

Requires `genesets-rs` built (`cargo build --release`). The GO+GOA download and
table build happen once per species and are cached under `--cache-dir`.

## What you do (judgment)

Read `curation/README.md` and the drafted associations, then per set:

- **skip** it if the screening report flags `MISMATCH` (or eyeball the drafted
  terms and confirm the membership doesn't support the ground truth);
- for each kept association, set `category` (core_process / core_component /
  supporting_process / nonspecific / false_association / …), `confidence`,
  `specificity`, and `insight` — the drafted `overlap_genes` are your grounded
  carrier evidence for the `curator_note`;
- **confirm or replace** the `AUTO-SUGGESTED` context term; add anatomy/cell-type
  contexts as needed;
- prune the long tail of nonspecific terms; add any `membership_gap` core terms
  the enrichment could not recover (curator-added);
- move the file to `curation/cgsa/genesets/`, then `curate validate` it.

The tool never sets `category`, `insight`, or a final context — those are the
curator's, by design. It removes the deterministic drudgery (membership lookup,
carrier grounding, marker/mismatch screening, OLS id verification via GOA), not
the interpretation.
