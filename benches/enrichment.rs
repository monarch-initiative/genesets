use std::hint::black_box;

use criterion::{Criterion, criterion_group, criterion_main};
use genesets_rs::{
    Correction, EnrichmentOptions,
    bitset::DenseBitSet,
    enrichment::enrich_matrix,
    index::{GeneUniverse, NamedSet, union_sets},
    output::{OutputFormat, OutputOptions, write_parquet_rows, write_rows},
};

fn synthetic_sets(
    count: usize,
    genes: usize,
    stride: usize,
    width: usize,
    prefix: &str,
) -> Vec<NamedSet> {
    (0..count)
        .map(|set_index| {
            let mut bits = DenseBitSet::new(genes);
            let start = (set_index * stride) % genes;
            for offset in 0..width {
                bits.set((start + offset * 17 + set_index) % genes);
            }
            NamedSet {
                id: format!("{prefix}_{set_index}"),
                name: None,
                bits,
            }
        })
        .collect()
}

fn bench_matrix(c: &mut Criterion) {
    let gene_count = 20_000;
    let targets = synthetic_sets(750, gene_count, 23, 300, "target");
    let queries = synthetic_sets(150, gene_count, 97, 250, "query");
    let background = union_sets(&targets, gene_count);
    let options = EnrichmentOptions {
        min_overlap: 1,
        correction: Correction::Bonferroni,
        include_overlap_genes: false,
        max_p_value: None,
        max_p_adjust: None,
    };

    c.bench_function("matrix_150x750_20k_genes", |b| {
        b.iter(|| black_box(enrich_matrix(&queries, &targets, &background, &options)))
    });

    let rows = enrich_matrix(&queries, &targets, &background, &options);
    let universe = GeneUniverse::new();
    c.bench_function("write_tsv_150x750_rows", |b| {
        b.iter(|| {
            let mut buffer = Vec::new();
            write_rows(
                &mut buffer,
                &rows,
                &queries,
                &targets,
                &universe,
                OutputOptions {
                    format: OutputFormat::Tsv,
                    include_overlap_genes: false,
                    correction: Correction::Bonferroni,
                },
            )
            .unwrap();
            black_box(buffer.len())
        })
    });

    c.bench_function("write_null_150x750_rows", |b| {
        b.iter(|| {
            let mut buffer = Vec::new();
            write_rows(
                &mut buffer,
                &rows,
                &queries,
                &targets,
                &universe,
                OutputOptions {
                    format: OutputFormat::Null,
                    include_overlap_genes: false,
                    correction: Correction::Bonferroni,
                },
            )
            .unwrap();
            black_box(buffer.len())
        })
    });

    c.bench_function("write_parquet_150x750_rows", |b| {
        b.iter(|| {
            let mut buffer = Vec::new();
            write_parquet_rows(
                &mut buffer,
                &rows,
                &queries,
                &targets,
                &universe,
                OutputOptions {
                    format: OutputFormat::Parquet,
                    include_overlap_genes: false,
                    correction: Correction::Bonferroni,
                },
            )
            .unwrap();
            black_box(buffer.len())
        })
    });
}

criterion_group!(benches, bench_matrix);
criterion_main!(benches);
