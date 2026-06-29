# Static Browser Plan

The curated interpretation corpus should become browsable from docs or a
dedicated static site before it needs a deployed dynamic service.

## Recommended Direction

Generate static pages from `curation/genesets/*.yaml`:

```text
curation/genesets/*.yaml
  -> validate with LinkML, term validator, and reference validator
  -> materialize a small JSON index
  -> generate one page per curated gene set
  -> generate collection, context, category, and recovery-status index pages
```

The generator should live in the Python workflow package, for example as:

```bash
genesets-workflows curate build-pages \
  --dir curation/genesets \
  --out-dir docs/src/generated/curation
```

Those generated pages can be committed while the corpus is small, or built in
CI for a static docs deploy once the corpus grows.

## Page Shape

Each gene set page should show:

- source ID, name, collection, taxon, direction, and gene count;
- biological contexts with ontology IDs and labels;
- curated GO associations grouped by role category;
- confidence, specificity, and recovery status;
- enrichment stats and overlap genes when present;
- evidence references and curator notes;
- links to neighboring index pages by collection, context, GO term, category,
  and recovery status.

## LinkML Browser Role

Use generated LinkML docs for the schema reference. Use generated corpus pages
for the curated instances. A schema browser explains what `TermAssociation`
means; an instance browser explains what
`MSIGDB:HAY_BONE_MARROW_DENDRITIC_CELL` means.

## Web UI Relationship

The current web explorer browses eval report bundles. A future curated gene set
browser can reuse the same principles but should read curated corpus artifacts
instead of report bundles:

- source of truth remains the validated YAML corpus;
- generated JSON gives the browser a stable read model;
- static pages support simple hosting and durable links;
- a richer client-side UI can be added later over the same JSON.

Do not add curated-browser behavior to the Rust compute engine. If enrichment
is needed to refresh stats, run that through workflows and write artifacts back
to the curation/report layer.
