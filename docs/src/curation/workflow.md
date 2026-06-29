# Curation Workflow

The curation workflow turns enrichment output into LinkML-validated biological
interpretations of non-GO gene sets.

## Local Setup

Use the curation extra from the Python workflow package:

```bash
uv run --project python/genesets-workflows --extra curation \
  genesets-workflows curate --help
```

The `justfile` provides the common repository commands:

```bash
just curate-validate-schema
just curate-validate
just curate-report
```

## Steps

1. Validate the schema's ontology-backed enum meanings:

   ```bash
   just curate-validate-schema
   ```

2. Draft a curated interpretation from enrichment output:

   ```bash
   uv run --project python/genesets-workflows --extra curation \
     genesets-workflows curate draft MSIGDB:<SET> \
     --enrichment-tsv path/to/enrichment.tsv \
     -o curation/genesets/<SET>.yaml
   ```

3. Adjudicate the draft by assigning `category`, `confidence`,
   `specificity`, and `recovery_status`. Add curator-only core terms and
   literature evidence where needed.

4. Validate every curated YAML file:

   ```bash
   just curate-validate
   ```

5. Build the precision/recall report:

   ```bash
   just curate-report
   ```

## Judgment Model

`category` is the biological judgment. A term that is central to the gene set's
biology should remain core even if current GO annotations do not recover it.

`recovery_status` records whether the current annotations and membership make
that biological judgment recoverable:

- `annotation_supported`: current annotations recover the term.
- `annotation_gap`: relevant genes are present, but GO annotation is too shallow.
- `membership_gap`: the gene set lacks genes needed to support the term.

This separation keeps the corpus useful both as an enrichment eval fixture and
as a source of curation targets.
