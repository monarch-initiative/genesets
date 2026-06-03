# Expression500 vs Reactome Flat

This eval enriches the same 500 expression-derived MSigDB/GSE query sets used
by the GO impact report against the official Reactome pathway GMT as a flat
target library.

Prepare the Reactome target library:

```bash
genesets-workflows prepare-reactome-flat
```

Run the enrichment:

```bash
genesets-rs run evals/reactome_flat/config.yaml
```

Or run both steps and write companion metadata:

```bash
genesets-workflows reactome-flat
```

The prep step normalizes Reactome's GMT from:

```text
pathway name    pathway stable id    gene symbols...
```

to the genesets-rs convention:

```text
pathway stable id    pathway name    gene symbols...
```

This is intentionally flat. A hierarchy-aware Reactome workflow can later use
the same core engine by preparing `terms.tsv`, `closure.tsv`, and
`gene_terms.tsv` from Reactome pathway relation and mapping files.
