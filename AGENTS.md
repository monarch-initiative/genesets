# AGENTS.md — orientation for agents and contributors

Guidance for working in this repo. Read this first, then the skills and READMEs it
points to. (Applies to any coding agent; Claude Code also auto-loads the skills in
`.claude/skills/`.)

## What this repo is

- **`genesets-rs`** (Rust, `src/`, `Cargo.toml`) — a fast Fisher-exact GO
  enrichment tool. Build with `cargo build --release`; the binary is `genesets-rs`.
- **`python/genesets-workflows`** — a Python (uv) workflow + curation layer with a
  `genesets-workflows` CLI (typer). Run anything in it via
  `uv run --project python/genesets-workflows --extra curation …`.
- **`curation/`** — a hand-curated **GO-interpretation gold standard**: one
  `genesets/<SET>.yaml` per non-GO gene set recording, per GO term, a curator
  judgment. It doubles as a precision/recall benchmark for the enrichment tool.
- **`evals/`** — evaluations, including `iba_vs_benchmark/` (scoring GOA evidence
  variants against the curated gold).

## The two workflows you'll most likely continue

| Task | Skill | Deep doc |
|---|---|---|
| Add/review a curated GO interpretation of a gene set | `.claude/skills/curate-geneset/SKILL.md` | `curation/README.md` + `curation/schema/genesets_interpretation.yaml` |
| Evaluate an enrichment method / GOA evidence variant vs the gold | `.claude/skills/evaluate-enrichment/SKILL.md` | `evals/iba_vs_benchmark/README.md` |

The skills are runbooks; the READMEs are the authority for the schema, the
categories, and the findings. Don't duplicate — extend.

## Conventions (this repo)

- **uv, never `requirements.txt`.** Python commands go through
  `uv run --project python/genesets-workflows …`.
- **`just` over raw commands.** Key targets: `just curate-validate` (4-gate
  validation of every interpretation), `just curate-test` (pytest + doctests),
  `just curate-report` (precision/recall/F1), `just curate-validate-schema`.
- **No try/except cheating.** If something throws unexpectedly, fix the cause;
  don't wrap it to swallow the error.
- **Doctests are tests** — keep them green (`just curate-test`).

## The non-negotiable curation discipline (why this gold standard is trustworthy)

- **Ground in the real membership.** Curate a set from its actual gene list, not
  the cell type's textbook markers. This is what makes `recovery_status` honest.
- **OLS-verify every ontology id+label and reject obsolete terms.** The validator's
  4th gate sweeps OAK `obsoletes()`; the `sqlite:obo` builds retain labels for
  obsolete classes, so id+label validation alone is not enough.
- **Evidence snippets are verbatim or absent.** The reference-validator
  substring-checks every `snippet` against the cited paper. Never fabricate a
  quote, a gene, or a PMID.
- **The eval measures the gold; it must never refit it.** `category` (biology) is
  the scored authority. `recovery_status`/`insight` are curator judgments
  adjudicated on GOA facts, never auto-updated to match a tool's output — that
  would make recall-vs-gold circular. See the Guardrail section in
  `evals/iba_vs_benchmark/README.md`. **If you are editing the gold to raise a
  recall number, stop.**

## Before you open / extend a PR

1. `just curate-validate` — all interpretation YAMLs pass the 4 gates.
2. `just curate-test` — unit + doctests pass.
3. Regenerate `curation/genesets/manifest.tsv` if you added/changed sets.
4. If you added evaluable sets, re-run the **evaluate-enrichment** pipeline and
   refresh the numbers in the eval README.

## Schema axes at a glance (full detail in `curation/README.md`)

- `category` — the biology (core/supporting/nonspecific/…); authoritative.
- `recovery_status` — `annotation_supported` / `annotation_gap` / `membership_gap`;
  grounded in the actual membership; orthogonal to `category`.
- `insight` — `confirmatory` vs `mechanistic`; the eval's recall splits on this
  (mechanistic insight is ~2× harder to recover).
- `series` / `series_role` — link contrasting poles of one axis (both poles get the
  field) so the eval can check they resolve to contrasting interpretations.

## Sources in the corpus

MSigDB (`MSIGDB:` — C8 single-cell, C2 pathways, H hallmark, C5:HPO phenotypes)
and literature (`LIT:` — `DISEASE_ACTIVITY` cell-state programs and `GENETIC`
GWAS/CRISPR convergences, with membership captured into
`curation/genesets/lit_members.gmt`).

## Gotchas

- `origin` (`cmungall/genesets-rs`) **redirects to `monarch-initiative/genesets`**;
  pushes print a "repository moved" notice but succeed, and PRs land there.
- The Descartes fetal atlas (`DESCARTES_FETAL_*` sets) is **PMID:33184181**, not
  `32848094` (a Dravet-syndrome paper) — a citation that's easy to get wrong.
