# Input Model

The engine has one central abstraction: a named set of genes. Ontology terms, flat database entries, user samples, and backgrounds can all be represented this way.

## Ontology Targets

Ontology enrichment targets are built from three tables:

Terms:

```text
term_id    name
GO:0006915 apoptosis
```

Closure:

```text
child      ancestor
GO:0006915 GO:0006915
GO:0006915 GO:0012501
```

Annotations:

```text
gene_id    term_id
TP53       GO:0006915
```

The closure table is expected to be reflexive, but the loader defensively adds each term as its own ancestor. If a gene is annotated to a child, the prepared bitset for every ancestor receives that gene.

## Flat Targets

Flat targets skip closure. They can be read from:

- `list`: one gene per line, one set total;
- `pairwise`: `set_id`, `gene_id`;
- `gene-term`: `gene_id`, `set_id`;
- `gmt`: `set_id`, `description`, genes;
- `gmx`: first row set IDs, remaining rows genes;
- `gmx-desc`: first row set IDs, second row descriptions, remaining rows genes.

## Background

If no background file is supplied, the background is the union of all target genes. A supplied background file overrides that behavior and is read as one gene per row.

This matters for experimental data. For RNA-seq differential-expression lists, the background should usually be all genes that were measured and eligible to become significant, not every protein-coding gene.

## Gene Names

Gene names are optional. The statistical engine uses gene IDs only. A gene-name table is only used for human-readable overlap output.
