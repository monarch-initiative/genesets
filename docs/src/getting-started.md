# Getting Started

Build the CLI:

```bash
cargo build --release
```

Install this checkout so `genesets-rs` can be run from any directory:

```bash
cargo install --path /path/to/genesets-rs --force
```

Cargo writes the binary to `~/.cargo/bin/genesets-rs`. If your shell cannot
find it, add Cargo's bin directory to your `PATH`:

```bash
echo 'export PATH="$HOME/.cargo/bin:$PATH"' >> ~/.zshrc
source ~/.zshrc
```

Run the bundled ontology-style example:

```bash
genesets-rs enrich \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --sample examples/sample.txt \
  --background examples/background.txt \
  --overlap-genes
```

Run all ontology terms against all ontology terms:

```bash
genesets-rs matrix \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --queries-from-targets \
  --background examples/background.txt
```

Run from YAML:

```bash
genesets-rs run examples/enrich.yaml
```

Relative paths inside a YAML config are resolved from the config file's
directory, so this also works from outside the repository:

```bash
genesets-rs run /path/to/genesets-rs/examples/enrich.yaml
```

For development without installing, use `cargo run --` before the same CLI
arguments:

```bash
cargo run -- run examples/enrich.yaml
```

Install or run the optional workflow layer when you want repeatable source
prep, GO/Reactome evals, Parquet summaries, and report metadata:

```bash
uv run --project python/genesets-workflows genesets-workflows doctor
```

The workflow layer calls the Rust CLI for batch compute; it is not required for
simple one-off enrichment.

Build the documentation site locally:

```bash
cargo install mdbook
mdbook serve
```

Run the test suite:

```bash
cargo test
```

Run the benchmark harness:

```bash
cargo bench
```
