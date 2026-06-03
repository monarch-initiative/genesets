# Storage Backends

The result backend should optimize for fast writes, compact artifacts, and easy
diff queries over many runs. TSV should remain available for small, inspectable
single analyses, but it should not be the mass-eval format.

## Recommendation

Use Parquet as the first durable mass-eval output format. Use DuckDB as the
default analysis layer over Parquet files rather than making `.duckdb` the first
primary output format.

This gives us:

- one fast, compact artifact format for CI, eval archives, and sharing;
- direct SQL introspection through DuckDB with no import step;
- simpler Rust output code using Arrow record batches and the native Parquet
  writer;
- a clean path to add a DuckDB sink later if repeated interactive sessions need
  a materialized database.

Supporting both as eventual output options is reasonable, but implementing both
as first-class write paths immediately is likely overkill. The first target
should be:

```text
genesets-rs ... --output-format parquet --output results.parquet
duckdb -c "SELECT * FROM 'results.parquet' WHERE p_adjust_bonferroni <= 0.05"
```

Diff outputs use the same approach:

```text
genesets-rs compare --left old.parquet --right new.parquet \
  --output-format parquet --output old-vs-new.diff.parquet
duckdb -c "SELECT class, count(*) FROM 'old-vs-new.diff.parquet' GROUP BY class"
```

## Comparative Notes

| Criterion | Parquet primary | DuckDB primary |
| --- | --- | --- |
| Rust write path | Direct Arrow `RecordBatch` to Parquet writer | `duckdb-rs` connection/appender or Arrow append |
| Artifact shape | Immutable columnar file or partitioned dataset | Embedded analytical database file |
| Queryability | Query directly with DuckDB, Polars, Arrow, Spark, Python | Query directly with DuckDB |
| File size | Usually smallest, especially with zstd | Slightly larger because it stores database/catalog structure |
| Append workflow | Prefer partitioned files, not in-place append | Natural append into tables |
| Reproducibility | Strong fit for immutable eval artifacts | Good, but easier to mutate accidentally |
| Distribution | Simple files/directories | One DB file, but DuckDB-version coupling matters more |
| Best use | Batch eval outputs and archived comparisons | Interactive workspaces and repeatedly queried derived tables |

DuckDB can read Parquet directly and can push projections and filters into the
Parquet scan. That means we do not need to choose between Parquet artifacts and
DuckDB introspection.

## Local Smoke Test

Using the 200 MSigDB GSE-style query fixture against current GOA human, the
unfiltered TSV had 778,450 result rows.

| Artifact | Size |
| --- | ---: |
| TSV | 132 MB |
| DuckDB table database | 21 MB |
| Parquet, DuckDB snappy export | 18 MB |
| Parquet, DuckDB zstd level 1 export | 11 MB |
| Parquet, `genesets-rs` snappy output | 11 MB |

DuckDB query timings, including process startup, were effectively tied on this
fixture:

| Query | DuckDB table | zstd Parquet read by DuckDB |
| --- | ---: | ---: |
| count significant rows | 0.03s | 0.03s |
| group significant rows by query | 0.04s | 0.03s |

These numbers are too small to settle large-scale behavior, but they are enough
to reject the idea that Parquet would cost us DuckDB-style introspection.

The first Rust Parquet backend is optimized for compatibility and write speed
with Snappy compression. On the same fixture, end-to-end unfiltered output was
about 3.7s for Parquet versus about 3.5s for TSV, while reducing output size
from 132 MB to 11 MB. The synthetic writer benchmark is more favorable to
Parquet: about 1.45 ms for Parquet versus about 3.7 ms for TSV on the current
150 x 750 fixture.

## Implementation Direction

The output boundary should produce numeric, batched result records:

- `run_id`, `query_index`, `target_index`;
- `overlap`, `query_size`, `target_size`, `background_size`;
- `p_value`, `p_adjust`, and boolean significance flags;
- optional string dictionaries or dimension tables for query and target labels;
- optional overlap genes in a separate sidecar table/file.

The Parquet backend writes Snappy-compressed row groups from these batches. A
later DuckDB backend can consume the same batches through an appender or
materialize a DuckDB database from the Parquet outputs.
