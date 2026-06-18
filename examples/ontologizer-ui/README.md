# Ontologizer UI Smoke Test Inputs

These files are plain-text gene lists for manually testing the P2GX Ontologizer UI:

- `population_genes_human.txt`: use as the Population Gene File.
- `study_genes_human.txt`: use as the Study Gene File.

They were copied from the MIT-licensed `P2GX/ontologizer` example data under
`examples/human/` and are intended as a small reproducible fixture for checking
that the GUI can load human GOA-style gene symbols and run both frequentist and
Bayesian analyses.

Use them with a human GAF file, such as `goa_human.gaf` or `goa_human.gaf.gz`,
and a matching GO ontology file, such as `go-basic.json`.

When run in Bayesian/MGSA mode, expected high-scoring terms include:

- `GO:0006513` protein monoubiquitination
- `GO:0006937` regulation of muscle contraction
- `GO:0030018` Z disc
- `GO:0031032` actomyosin structure organization
- `GO:0031648` protein destabilization
- `GO:0032663` regulation of interleukin-2 production
- `GO:0045095` keratin filament
- `GO:0072332` intrinsic apoptotic signaling pathway by p53 class mediator
- `GO:0140359` ABC-type transporter activity

Exact labels, study-hit lists, and posterior probabilities can drift with the
GO and GOA release used by the UI.
