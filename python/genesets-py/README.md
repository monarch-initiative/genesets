# genesets (Python)

Python bindings for [genesets-rs](https://github.com/monarch-initiative/genesets).
The compute runs in Rust and releases the GIL.

```python
import genesets

result = genesets.mgsa(
    annotations=[("BRCA1", "GO:0006281"), ...],   # (gene, term), propagated to ancestors
    study=["BRCA1", "RAD51", ...],
    background=None,                               # default: every annotated gene
    closure=[("GO:0006281", "GO:0006974"), ...],   # optional (child, ancestor) pairs
)
for term in result["terms"][:10]:
    print(term["id"], term["posterior"], term["overlap"], term["term_size"])
```

`mgsa` is Bayesian model-based gene set analysis, as in Ontologizer and
wtfgenes. Each term gets a posterior probability of being active, so ancestors
of an enriched term are explained away instead of being reported alongside it.
Keyword options match the `genesets-rs mgsa` CLI (`samples_per_term`,
`burn_per_term`, the `*_prob`/`*_count` priors, move rates, and `seed`).

## Install

Until wheels are on PyPI, install from git. Building from source needs a Rust
toolchain:

```bash
uv add "genesets @ git+https://github.com/monarch-initiative/genesets#subdirectory=python/genesets-py"
```

## Develop

```bash
cd python/genesets-py
uv run --with maturin --with pytest maturin develop --release
uv run --with pytest pytest tests
```
