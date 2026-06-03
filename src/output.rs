use std::{io::Write, sync::Arc};

use anyhow::{Context, Result, bail};
use arrow_array::{
    ArrayRef, BooleanArray, Float64Array, Int32Array, RecordBatch, StringArray, UInt32Array,
};
use arrow_schema::{DataType, Field, Schema, SchemaRef};
use clap::ValueEnum;
use parquet::{arrow::ArrowWriter, basic::Compression, file::properties::WriterProperties};
use serde::Deserialize;

use crate::{
    compare::DiffRow,
    enrichment::{Correction, EnrichmentRow},
    index::{GeneUniverse, NamedSet},
};

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[value(rename_all = "kebab-case")]
pub enum OutputFormat {
    Tsv,
    Parquet,
    Null,
}

#[derive(Copy, Clone, Debug)]
pub struct OutputOptions {
    pub format: OutputFormat,
    pub include_overlap_genes: bool,
    pub correction: Correction,
}

pub fn write_rows<W: Write + ?Sized>(
    writer: &mut W,
    rows: &[EnrichmentRow],
    queries: &[NamedSet],
    targets: &[NamedSet],
    genes: &GeneUniverse,
    options: OutputOptions,
) -> Result<()> {
    match options.format {
        OutputFormat::Tsv => write_tsv_rows(
            writer,
            rows,
            queries,
            targets,
            genes,
            options.include_overlap_genes,
            options.correction,
        ),
        OutputFormat::Parquet => {
            bail!("parquet output requires a file path; use --output with --output-format parquet")
        }
        OutputFormat::Null => Ok(()),
    }
}

pub fn write_parquet_rows<W: Write + Send>(
    writer: W,
    rows: &[EnrichmentRow],
    queries: &[NamedSet],
    targets: &[NamedSet],
    genes: &GeneUniverse,
    options: OutputOptions,
) -> Result<()> {
    let schema = parquet_schema(options.include_overlap_genes, options.correction);
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_row_count(Some(65_536))
        .set_created_by(format!("genesets-rs {}", env!("CARGO_PKG_VERSION")))
        .build();
    let mut writer = ArrowWriter::try_new(writer, Arc::clone(&schema), Some(props))?;

    for chunk in rows.chunks(65_536) {
        let batch = build_record_batch(chunk, queries, targets, genes, &schema, options)?;
        writer.write(&batch)?;
    }

    writer.close()?;
    Ok(())
}

pub fn write_diff_rows<W: Write + ?Sized>(writer: &mut W, rows: &[DiffRow]) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    writer.write_record([
        "class",
        "query_id",
        "target_id",
        "query_name",
        "target_name",
        "left_present",
        "right_present",
        "left_significant",
        "right_significant",
        "left_p_adjust",
        "right_p_adjust",
        "delta_neg_log10_p_adjust",
        "left_p_value",
        "right_p_value",
        "left_overlap",
        "right_overlap",
        "delta_overlap",
        "left_query_size",
        "right_query_size",
        "left_target_size",
        "right_target_size",
        "delta_target_size",
        "left_background_size",
        "right_background_size",
    ])?;
    for row in rows {
        writer.write_record([
            row.class.as_str().to_owned(),
            row.query_id.clone(),
            row.target_id.clone(),
            row.query_name.clone(),
            row.target_name.clone(),
            row.left_present.to_string(),
            row.right_present.to_string(),
            row.left_significant.to_string(),
            row.right_significant.to_string(),
            format_optional_f64(row.left_p_adjust),
            format_optional_f64(row.right_p_adjust),
            format_optional_f64(row.delta_neg_log10_p_adjust),
            format_optional_f64(row.left_p_value),
            format_optional_f64(row.right_p_value),
            format_optional_u32(row.left_overlap),
            format_optional_u32(row.right_overlap),
            format_optional_i32(row.delta_overlap),
            format_optional_u32(row.left_query_size),
            format_optional_u32(row.right_query_size),
            format_optional_u32(row.left_target_size),
            format_optional_u32(row.right_target_size),
            format_optional_i32(row.delta_target_size),
            format_optional_u32(row.left_background_size),
            format_optional_u32(row.right_background_size),
        ])?;
    }
    writer.flush()?;
    Ok(())
}

pub fn write_parquet_diff_rows<W: Write + Send>(writer: W, rows: &[DiffRow]) -> Result<()> {
    let schema = diff_parquet_schema();
    let props = WriterProperties::builder()
        .set_compression(Compression::SNAPPY)
        .set_max_row_group_row_count(Some(65_536))
        .set_created_by(format!("genesets-rs {}", env!("CARGO_PKG_VERSION")))
        .build();
    let mut writer = ArrowWriter::try_new(writer, Arc::clone(&schema), Some(props))?;
    for chunk in rows.chunks(65_536) {
        let batch = build_diff_record_batch(chunk, &schema)?;
        writer.write(&batch)?;
    }
    writer.close()?;
    Ok(())
}

fn write_tsv_rows<W: Write + ?Sized>(
    writer: &mut W,
    rows: &[EnrichmentRow],
    queries: &[NamedSet],
    targets: &[NamedSet],
    genes: &GeneUniverse,
    include_overlap_genes: bool,
    correction: Correction,
) -> Result<()> {
    let mut writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    let correction_column = match correction {
        Correction::Bonferroni => "p_adjust_bonferroni",
        Correction::None => "p_adjust",
    };

    let mut header = vec![
        "query_id",
        "query_name",
        "target_id",
        "target_name",
        "overlap",
        "query_size",
        "target_size",
        "background_size",
        "p_value",
        correction_column,
    ];
    if include_overlap_genes {
        header.push("overlap_genes");
        header.push("overlap_gene_names");
    }
    writer.write_record(header)?;

    for row in rows {
        let query = &queries[row.query_index];
        let target = &targets[row.target_index];
        let mut record = vec![
            query.id.clone(),
            query.name.clone().unwrap_or_default(),
            target.id.clone(),
            target.name.clone().unwrap_or_default(),
            row.overlap.to_string(),
            row.query_size.to_string(),
            row.target_size.to_string(),
            row.background_size.to_string(),
            format_p_value(row.p_value),
            format_p_value(row.p_adjust),
        ];

        if include_overlap_genes {
            let indices = row.overlap_indices.as_deref().unwrap_or(&[]);
            let gene_ids = indices
                .iter()
                .map(|index| genes.gene_id(*index))
                .collect::<Vec<_>>()
                .join(";");
            let gene_names = indices
                .iter()
                .filter_map(|index| genes.gene_name(*index))
                .collect::<Vec<_>>()
                .join(";");
            record.push(gene_ids);
            record.push(gene_names);
        }

        writer.write_record(record)?;
    }

    writer.flush()?;
    Ok(())
}

fn format_p_value(value: f64) -> String {
    if value == 0.0 {
        "0".to_owned()
    } else if value.is_nan() {
        "NaN".to_owned()
    } else {
        format!("{value:.6e}")
    }
}

fn format_optional_f64(value: Option<f64>) -> String {
    value.map(format_p_value).unwrap_or_default()
}

fn format_optional_u32(value: Option<u32>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn format_optional_i32(value: Option<i32>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn parquet_schema(include_overlap_genes: bool, correction: Correction) -> SchemaRef {
    let correction_column = match correction {
        Correction::Bonferroni => "p_adjust_bonferroni",
        Correction::None => "p_adjust",
    };
    let mut fields = vec![
        Field::new("query_index", DataType::UInt32, false),
        Field::new("target_index", DataType::UInt32, false),
        Field::new("query_id", DataType::Utf8, false),
        Field::new("query_name", DataType::Utf8, false),
        Field::new("target_id", DataType::Utf8, false),
        Field::new("target_name", DataType::Utf8, false),
        Field::new("overlap", DataType::UInt32, false),
        Field::new("query_size", DataType::UInt32, false),
        Field::new("target_size", DataType::UInt32, false),
        Field::new("background_size", DataType::UInt32, false),
        Field::new("p_value", DataType::Float64, false),
        Field::new(correction_column, DataType::Float64, false),
    ];
    if include_overlap_genes {
        fields.push(Field::new("overlap_genes", DataType::Utf8, false));
        fields.push(Field::new("overlap_gene_names", DataType::Utf8, false));
    }
    Arc::new(Schema::new(fields))
}

fn build_record_batch(
    rows: &[EnrichmentRow],
    queries: &[NamedSet],
    targets: &[NamedSet],
    genes: &GeneUniverse,
    schema: &SchemaRef,
    options: OutputOptions,
) -> Result<RecordBatch> {
    let mut query_indices = Vec::with_capacity(rows.len());
    let mut target_indices = Vec::with_capacity(rows.len());
    let mut query_ids = Vec::with_capacity(rows.len());
    let mut query_names = Vec::with_capacity(rows.len());
    let mut target_ids = Vec::with_capacity(rows.len());
    let mut target_names = Vec::with_capacity(rows.len());
    let mut overlaps = Vec::with_capacity(rows.len());
    let mut query_sizes = Vec::with_capacity(rows.len());
    let mut target_sizes = Vec::with_capacity(rows.len());
    let mut background_sizes = Vec::with_capacity(rows.len());
    let mut p_values = Vec::with_capacity(rows.len());
    let mut p_adjusts = Vec::with_capacity(rows.len());
    let mut overlap_genes = Vec::with_capacity(rows.len());
    let mut overlap_gene_names = Vec::with_capacity(rows.len());

    for row in rows {
        let query = &queries[row.query_index];
        let target = &targets[row.target_index];
        query_indices
            .push(u32::try_from(row.query_index).context("query index does not fit in UInt32")?);
        target_indices
            .push(u32::try_from(row.target_index).context("target index does not fit in UInt32")?);
        query_ids.push(query.id.clone());
        query_names.push(query.name.clone().unwrap_or_default());
        target_ids.push(target.id.clone());
        target_names.push(target.name.clone().unwrap_or_default());
        overlaps.push(row.overlap);
        query_sizes.push(row.query_size);
        target_sizes.push(row.target_size);
        background_sizes.push(row.background_size);
        p_values.push(row.p_value);
        p_adjusts.push(row.p_adjust);

        if options.include_overlap_genes {
            let indices = row.overlap_indices.as_deref().unwrap_or(&[]);
            overlap_genes.push(
                indices
                    .iter()
                    .map(|index| genes.gene_id(*index))
                    .collect::<Vec<_>>()
                    .join(";"),
            );
            overlap_gene_names.push(
                indices
                    .iter()
                    .filter_map(|index| genes.gene_name(*index))
                    .collect::<Vec<_>>()
                    .join(";"),
            );
        }
    }

    let mut columns: Vec<ArrayRef> = vec![
        Arc::new(UInt32Array::from(query_indices)),
        Arc::new(UInt32Array::from(target_indices)),
        Arc::new(StringArray::from(query_ids)),
        Arc::new(StringArray::from(query_names)),
        Arc::new(StringArray::from(target_ids)),
        Arc::new(StringArray::from(target_names)),
        Arc::new(UInt32Array::from(overlaps)),
        Arc::new(UInt32Array::from(query_sizes)),
        Arc::new(UInt32Array::from(target_sizes)),
        Arc::new(UInt32Array::from(background_sizes)),
        Arc::new(Float64Array::from(p_values)),
        Arc::new(Float64Array::from(p_adjusts)),
    ];
    if options.include_overlap_genes {
        columns.push(Arc::new(StringArray::from(overlap_genes)));
        columns.push(Arc::new(StringArray::from(overlap_gene_names)));
    }

    Ok(RecordBatch::try_new(Arc::clone(schema), columns)?)
}

fn diff_parquet_schema() -> SchemaRef {
    Arc::new(Schema::new(vec![
        Field::new("class", DataType::Utf8, false),
        Field::new("query_id", DataType::Utf8, false),
        Field::new("target_id", DataType::Utf8, false),
        Field::new("query_name", DataType::Utf8, false),
        Field::new("target_name", DataType::Utf8, false),
        Field::new("left_present", DataType::Boolean, false),
        Field::new("right_present", DataType::Boolean, false),
        Field::new("left_significant", DataType::Boolean, false),
        Field::new("right_significant", DataType::Boolean, false),
        Field::new("left_p_adjust", DataType::Float64, true),
        Field::new("right_p_adjust", DataType::Float64, true),
        Field::new("delta_neg_log10_p_adjust", DataType::Float64, true),
        Field::new("left_p_value", DataType::Float64, true),
        Field::new("right_p_value", DataType::Float64, true),
        Field::new("left_overlap", DataType::UInt32, true),
        Field::new("right_overlap", DataType::UInt32, true),
        Field::new("delta_overlap", DataType::Int32, true),
        Field::new("left_query_size", DataType::UInt32, true),
        Field::new("right_query_size", DataType::UInt32, true),
        Field::new("left_target_size", DataType::UInt32, true),
        Field::new("right_target_size", DataType::UInt32, true),
        Field::new("delta_target_size", DataType::Int32, true),
        Field::new("left_background_size", DataType::UInt32, true),
        Field::new("right_background_size", DataType::UInt32, true),
    ]))
}

fn build_diff_record_batch(rows: &[DiffRow], schema: &SchemaRef) -> Result<RecordBatch> {
    Ok(RecordBatch::try_new(
        Arc::clone(schema),
        vec![
            Arc::new(StringArray::from(
                rows.iter()
                    .map(|row| row.class.as_str())
                    .collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(
                rows.iter()
                    .map(|row| row.query_id.clone())
                    .collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(
                rows.iter()
                    .map(|row| row.target_id.clone())
                    .collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(
                rows.iter()
                    .map(|row| row.query_name.clone())
                    .collect::<Vec<_>>(),
            )),
            Arc::new(StringArray::from(
                rows.iter()
                    .map(|row| row.target_name.clone())
                    .collect::<Vec<_>>(),
            )),
            Arc::new(BooleanArray::from(
                rows.iter().map(|row| row.left_present).collect::<Vec<_>>(),
            )),
            Arc::new(BooleanArray::from(
                rows.iter().map(|row| row.right_present).collect::<Vec<_>>(),
            )),
            Arc::new(BooleanArray::from(
                rows.iter()
                    .map(|row| row.left_significant)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(BooleanArray::from(
                rows.iter()
                    .map(|row| row.right_significant)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(
                rows.iter().map(|row| row.left_p_adjust).collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(
                rows.iter()
                    .map(|row| row.right_p_adjust)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(
                rows.iter()
                    .map(|row| row.delta_neg_log10_p_adjust)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(
                rows.iter().map(|row| row.left_p_value).collect::<Vec<_>>(),
            )),
            Arc::new(Float64Array::from(
                rows.iter().map(|row| row.right_p_value).collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter().map(|row| row.left_overlap).collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter().map(|row| row.right_overlap).collect::<Vec<_>>(),
            )),
            Arc::new(Int32Array::from(
                rows.iter().map(|row| row.delta_overlap).collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.left_query_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.right_query_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.left_target_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.right_target_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(Int32Array::from(
                rows.iter()
                    .map(|row| row.delta_target_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.left_background_size)
                    .collect::<Vec<_>>(),
            )),
            Arc::new(UInt32Array::from(
                rows.iter()
                    .map(|row| row.right_background_size)
                    .collect::<Vec<_>>(),
            )),
        ],
    )?)
}
