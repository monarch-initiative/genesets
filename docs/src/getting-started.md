# Getting Started

Build the CLI:

```bash
cargo build --release
```

Run the bundled ontology-style example:

```bash
cargo run -- enrich \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --sample examples/sample.txt \
  --background examples/background.txt \
  --overlap-genes
```

Run all ontology terms against all ontology terms:

```bash
cargo run -- matrix \
  --annotations examples/gene_terms.tsv \
  --terms examples/terms.tsv \
  --closure examples/closure.tsv \
  --queries-from-targets \
  --background examples/background.txt
```

Run from YAML:

```bash
cargo run -- run examples/enrich.yaml
```

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
