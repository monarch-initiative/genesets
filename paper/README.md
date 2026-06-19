# Paper Draft

This directory contains a LaTeX manuscript draft for the `genesets-rs`
framework paper.

Build locally:

```bash
cd paper
latexmk -pdf main.tex
```

Clean build artifacts:

```bash
latexmk -C main.tex
```

The manuscript cites benchmark summaries generated under `notebooks/generated/`.
Large generated Parquet and GMT artifacts are intentionally not tracked in git.
