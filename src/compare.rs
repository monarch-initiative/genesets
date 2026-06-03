use std::{
    collections::HashMap,
    fs::File,
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use arrow_array::{
    Array, Float64Array, Int64Array, RecordBatch, StringArray, UInt32Array, UInt64Array,
};
use clap::ValueEnum;
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use serde::Serialize;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
#[value(rename_all = "kebab-case")]
pub enum ResultFormat {
    Auto,
    Tsv,
    Parquet,
}

#[derive(Clone, Debug)]
pub struct ResultRecord {
    pub query_id: String,
    pub query_name: String,
    pub target_id: String,
    pub target_name: String,
    pub overlap: Option<u32>,
    pub query_size: Option<u32>,
    pub target_size: Option<u32>,
    pub background_size: Option<u32>,
    pub p_value: Option<f64>,
    pub p_adjust: Option<f64>,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash)]
pub struct ResultKey {
    pub query_id: String,
    pub target_id: String,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub enum DiffClass {
    SharedSignificant,
    LostSignificant,
    GainedSignificant,
}

impl DiffClass {
    pub fn as_str(&self) -> &'static str {
        match self {
            Self::SharedSignificant => "shared_significant",
            Self::LostSignificant => "lost_significant",
            Self::GainedSignificant => "gained_significant",
        }
    }
}

#[derive(Clone, Debug)]
pub struct DiffOptions {
    pub p_adjust_cutoff: f64,
    pub crossings_only: bool,
}

#[derive(Clone, Debug)]
pub struct DiffRow {
    pub class: DiffClass,
    pub query_id: String,
    pub target_id: String,
    pub query_name: String,
    pub target_name: String,
    pub left_present: bool,
    pub right_present: bool,
    pub left_significant: bool,
    pub right_significant: bool,
    pub left_p_adjust: Option<f64>,
    pub right_p_adjust: Option<f64>,
    pub delta_neg_log10_p_adjust: Option<f64>,
    pub left_p_value: Option<f64>,
    pub right_p_value: Option<f64>,
    pub left_overlap: Option<u32>,
    pub right_overlap: Option<u32>,
    pub delta_overlap: Option<i32>,
    pub left_query_size: Option<u32>,
    pub right_query_size: Option<u32>,
    pub left_target_size: Option<u32>,
    pub right_target_size: Option<u32>,
    pub delta_target_size: Option<i32>,
    pub left_background_size: Option<u32>,
    pub right_background_size: Option<u32>,
}

#[derive(Clone, Debug, Serialize)]
pub struct DiffSummary {
    pub left_path: PathBuf,
    pub right_path: PathBuf,
    pub left_format: String,
    pub right_format: String,
    pub p_adjust_cutoff: f64,
    pub crossings_only: bool,
    pub left_rows: usize,
    pub right_rows: usize,
    pub output_rows: usize,
    pub shared_significant: usize,
    pub lost_significant: usize,
    pub gained_significant: usize,
}

#[derive(Clone, Debug)]
pub struct DiffSummaryContext {
    pub left_path: PathBuf,
    pub right_path: PathBuf,
    pub left_format: ResultFormat,
    pub right_format: ResultFormat,
    pub left_rows: usize,
    pub right_rows: usize,
}

pub fn load_result_records(path: &Path, format: ResultFormat) -> Result<Vec<ResultRecord>> {
    match resolve_format(path, format)? {
        ResultFormat::Tsv => read_tsv_result_records(path),
        ResultFormat::Parquet => read_parquet_result_records(path),
        ResultFormat::Auto => unreachable!("auto result format is resolved before loading"),
    }
}

pub fn compare_records(
    left: &[ResultRecord],
    right: &[ResultRecord],
    options: &DiffOptions,
) -> Vec<DiffRow> {
    let left_index = build_index(left);
    let right_index = build_index(right);
    let mut keys: Vec<&ResultKey> = left_index.keys().chain(right_index.keys()).collect();
    keys.sort_by(|left, right| {
        left.query_id
            .cmp(&right.query_id)
            .then_with(|| left.target_id.cmp(&right.target_id))
    });
    keys.dedup();

    let mut rows = Vec::new();
    for key in keys {
        let left = left_index.get(key).copied();
        let right = right_index.get(key).copied();
        let left_significant = is_significant(left, options.p_adjust_cutoff);
        let right_significant = is_significant(right, options.p_adjust_cutoff);
        let class = match (left_significant, right_significant) {
            (true, true) if !options.crossings_only => DiffClass::SharedSignificant,
            (true, false) => DiffClass::LostSignificant,
            (false, true) => DiffClass::GainedSignificant,
            _ => continue,
        };
        rows.push(build_diff_row(key, class, left, right));
    }

    rows.sort_by(|left, right| {
        class_rank(&left.class)
            .cmp(&class_rank(&right.class))
            .then_with(|| left.query_id.cmp(&right.query_id))
            .then_with(|| left.target_id.cmp(&right.target_id))
    });
    rows
}

pub fn summarize_diff(
    context: DiffSummaryContext,
    options: &DiffOptions,
    rows: &[DiffRow],
) -> DiffSummary {
    DiffSummary {
        left_path: context.left_path,
        right_path: context.right_path,
        left_format: format_name(context.left_format).to_owned(),
        right_format: format_name(context.right_format).to_owned(),
        p_adjust_cutoff: options.p_adjust_cutoff,
        crossings_only: options.crossings_only,
        left_rows: context.left_rows,
        right_rows: context.right_rows,
        output_rows: rows.len(),
        shared_significant: rows
            .iter()
            .filter(|row| row.class == DiffClass::SharedSignificant)
            .count(),
        lost_significant: rows
            .iter()
            .filter(|row| row.class == DiffClass::LostSignificant)
            .count(),
        gained_significant: rows
            .iter()
            .filter(|row| row.class == DiffClass::GainedSignificant)
            .count(),
    }
}

pub fn resolve_format(path: &Path, format: ResultFormat) -> Result<ResultFormat> {
    if format != ResultFormat::Auto {
        return Ok(format);
    }
    match path.extension().and_then(|extension| extension.to_str()) {
        Some("parquet" | "pq") => Ok(ResultFormat::Parquet),
        Some("tsv" | "txt") => Ok(ResultFormat::Tsv),
        _ => bail!(
            "could not infer result format for {}; pass --left-format or --right-format",
            path.display()
        ),
    }
}

pub fn format_name(format: ResultFormat) -> &'static str {
    match format {
        ResultFormat::Auto => "auto",
        ResultFormat::Tsv => "tsv",
        ResultFormat::Parquet => "parquet",
    }
}

fn read_tsv_result_records(path: &Path) -> Result<Vec<ResultRecord>> {
    let mut reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_path(path)
        .with_context(|| format!("failed to open result TSV {}", path.display()))?;
    let headers = reader
        .headers()
        .with_context(|| format!("failed to read TSV header from {}", path.display()))?
        .clone();
    let columns = ResultColumns::from_headers(headers.iter())?;
    let mut records = Vec::new();
    for row in reader.records() {
        let row = row.with_context(|| format!("failed to read row from {}", path.display()))?;
        records.push(ResultRecord {
            query_id: required_tsv_value(&row, columns.query_id, "query_id")?.to_owned(),
            query_name: optional_tsv_string(&row, columns.query_name),
            target_id: required_tsv_value(&row, columns.target_id, "target_id")?.to_owned(),
            target_name: optional_tsv_string(&row, columns.target_name),
            overlap: optional_tsv_u32(&row, columns.overlap)?,
            query_size: optional_tsv_u32(&row, columns.query_size)?,
            target_size: optional_tsv_u32(&row, columns.target_size)?,
            background_size: optional_tsv_u32(&row, columns.background_size)?,
            p_value: optional_tsv_f64(&row, columns.p_value)?,
            p_adjust: optional_tsv_f64(&row, Some(columns.p_adjust))?,
        });
    }
    Ok(records)
}

fn read_parquet_result_records(path: &Path) -> Result<Vec<ResultRecord>> {
    let file = File::open(path)
        .with_context(|| format!("failed to open result parquet {}", path.display()))?;
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)?
        .with_batch_size(65_536)
        .build()?;
    let mut records = Vec::new();
    for batch in reader {
        let batch =
            batch.with_context(|| format!("failed to read batch from {}", path.display()))?;
        let columns = ResultColumns::from_schema(&batch)?;
        for row_index in 0..batch.num_rows() {
            records.push(ResultRecord {
                query_id: string_value(&batch, columns.query_id, row_index, "query_id")?,
                query_name: optional_string_value(&batch, columns.query_name, row_index)?,
                target_id: string_value(&batch, columns.target_id, row_index, "target_id")?,
                target_name: optional_string_value(&batch, columns.target_name, row_index)?,
                overlap: optional_u32_value(&batch, columns.overlap, row_index)?,
                query_size: optional_u32_value(&batch, columns.query_size, row_index)?,
                target_size: optional_u32_value(&batch, columns.target_size, row_index)?,
                background_size: optional_u32_value(&batch, columns.background_size, row_index)?,
                p_value: optional_f64_value(&batch, columns.p_value, row_index)?,
                p_adjust: optional_f64_value(&batch, Some(columns.p_adjust), row_index)?,
            });
        }
    }
    Ok(records)
}

fn build_index(records: &[ResultRecord]) -> HashMap<ResultKey, &ResultRecord> {
    records
        .iter()
        .map(|record| {
            (
                ResultKey {
                    query_id: record.query_id.clone(),
                    target_id: record.target_id.clone(),
                },
                record,
            )
        })
        .collect()
}

fn is_significant(record: Option<&ResultRecord>, cutoff: f64) -> bool {
    record
        .and_then(|record| record.p_adjust)
        .is_some_and(|p_adjust| p_adjust <= cutoff)
}

fn build_diff_row(
    key: &ResultKey,
    class: DiffClass,
    left: Option<&ResultRecord>,
    right: Option<&ResultRecord>,
) -> DiffRow {
    let query_name = prefer_string(
        right.map(|record| record.query_name.as_str()),
        left.map(|record| record.query_name.as_str()),
    );
    let target_name = prefer_string(
        right.map(|record| record.target_name.as_str()),
        left.map(|record| record.target_name.as_str()),
    );
    let (left_significant, right_significant) = match class {
        DiffClass::SharedSignificant => (true, true),
        DiffClass::LostSignificant => (true, false),
        DiffClass::GainedSignificant => (false, true),
    };

    DiffRow {
        class,
        query_id: key.query_id.clone(),
        target_id: key.target_id.clone(),
        query_name,
        target_name,
        left_present: left.is_some(),
        right_present: right.is_some(),
        left_significant,
        right_significant,
        left_p_adjust: left.and_then(|record| record.p_adjust),
        right_p_adjust: right.and_then(|record| record.p_adjust),
        delta_neg_log10_p_adjust: match (
            left.and_then(|record| record.p_adjust),
            right.and_then(|record| record.p_adjust),
        ) {
            (Some(left), Some(right)) => Some(neg_log10(right) - neg_log10(left)),
            _ => None,
        },
        left_p_value: left.and_then(|record| record.p_value),
        right_p_value: right.and_then(|record| record.p_value),
        left_overlap: left.and_then(|record| record.overlap),
        right_overlap: right.and_then(|record| record.overlap),
        delta_overlap: diff_u32(
            left.and_then(|record| record.overlap),
            right.and_then(|record| record.overlap),
        ),
        left_query_size: left.and_then(|record| record.query_size),
        right_query_size: right.and_then(|record| record.query_size),
        left_target_size: left.and_then(|record| record.target_size),
        right_target_size: right.and_then(|record| record.target_size),
        delta_target_size: diff_u32(
            left.and_then(|record| record.target_size),
            right.and_then(|record| record.target_size),
        ),
        left_background_size: left.and_then(|record| record.background_size),
        right_background_size: right.and_then(|record| record.background_size),
    }
}

fn prefer_string(right: Option<&str>, left: Option<&str>) -> String {
    right
        .filter(|value| !value.is_empty())
        .or_else(|| left.filter(|value| !value.is_empty()))
        .unwrap_or_default()
        .to_owned()
}

fn diff_u32(left: Option<u32>, right: Option<u32>) -> Option<i32> {
    match (left, right) {
        (Some(left), Some(right)) => i32::try_from(i64::from(right) - i64::from(left)).ok(),
        _ => None,
    }
}

fn neg_log10(value: f64) -> f64 {
    if value <= 0.0 { 320.0 } else { -value.log10() }
}

fn class_rank(class: &DiffClass) -> u8 {
    match class {
        DiffClass::LostSignificant => 0,
        DiffClass::GainedSignificant => 1,
        DiffClass::SharedSignificant => 2,
    }
}

#[derive(Clone, Debug)]
struct ResultColumns {
    query_id: usize,
    query_name: Option<usize>,
    target_id: usize,
    target_name: Option<usize>,
    overlap: Option<usize>,
    query_size: Option<usize>,
    target_size: Option<usize>,
    background_size: Option<usize>,
    p_value: Option<usize>,
    p_adjust: usize,
}

impl ResultColumns {
    fn from_headers<'a>(headers: impl Iterator<Item = &'a str>) -> Result<Self> {
        let offsets: HashMap<&str, usize> =
            headers.enumerate().map(|(i, name)| (name, i)).collect();
        Self::from_offsets(&offsets)
    }

    fn from_schema(batch: &RecordBatch) -> Result<Self> {
        let schema = batch.schema();
        let offsets: HashMap<&str, usize> = schema
            .fields()
            .iter()
            .enumerate()
            .map(|(i, field)| (field.name().as_str(), i))
            .collect();
        Self::from_offsets(&offsets)
    }

    fn from_offsets(offsets: &HashMap<&str, usize>) -> Result<Self> {
        let p_adjust = offsets
            .get("p_adjust_bonferroni")
            .or_else(|| offsets.get("p_adjust"))
            .copied()
            .context("result table requires p_adjust_bonferroni or p_adjust")?;
        Ok(Self {
            query_id: required_offset(offsets, "query_id")?,
            query_name: offsets.get("query_name").copied(),
            target_id: required_offset(offsets, "target_id")?,
            target_name: offsets.get("target_name").copied(),
            overlap: offsets.get("overlap").copied(),
            query_size: offsets.get("query_size").copied(),
            target_size: offsets.get("target_size").copied(),
            background_size: offsets.get("background_size").copied(),
            p_value: offsets.get("p_value").copied(),
            p_adjust,
        })
    }
}

fn required_offset(offsets: &HashMap<&str, usize>, name: &str) -> Result<usize> {
    offsets
        .get(name)
        .copied()
        .with_context(|| format!("result table requires {name}"))
}

fn required_tsv_value<'a>(
    row: &'a csv::StringRecord,
    offset: usize,
    column: &str,
) -> Result<&'a str> {
    row.get(offset)
        .filter(|value| !value.is_empty())
        .with_context(|| format!("row is missing required column {column}"))
}

fn optional_tsv_string(row: &csv::StringRecord, offset: Option<usize>) -> String {
    offset
        .and_then(|offset| row.get(offset))
        .unwrap_or_default()
        .to_owned()
}

fn optional_tsv_f64(row: &csv::StringRecord, offset: Option<usize>) -> Result<Option<f64>> {
    let Some(value) = offset.and_then(|offset| row.get(offset)) else {
        return Ok(None);
    };
    if value.is_empty() {
        return Ok(None);
    }
    Ok(Some(value.parse()?))
}

fn optional_tsv_u32(row: &csv::StringRecord, offset: Option<usize>) -> Result<Option<u32>> {
    let Some(value) = offset.and_then(|offset| row.get(offset)) else {
        return Ok(None);
    };
    if value.is_empty() {
        return Ok(None);
    }
    Ok(Some(value.parse()?))
}

fn string_value(batch: &RecordBatch, offset: usize, row: usize, column: &str) -> Result<String> {
    let array = batch
        .column(offset)
        .as_any()
        .downcast_ref::<StringArray>()
        .with_context(|| format!("{column} must be Utf8"))?;
    if array.is_null(row) {
        bail!("{column} is null");
    }
    Ok(array.value(row).to_owned())
}

fn optional_string_value(batch: &RecordBatch, offset: Option<usize>, row: usize) -> Result<String> {
    let Some(offset) = offset else {
        return Ok(String::new());
    };
    let array = batch
        .column(offset)
        .as_any()
        .downcast_ref::<StringArray>()
        .context("optional string column must be Utf8")?;
    if array.is_null(row) {
        Ok(String::new())
    } else {
        Ok(array.value(row).to_owned())
    }
}

fn optional_f64_value(
    batch: &RecordBatch,
    offset: Option<usize>,
    row: usize,
) -> Result<Option<f64>> {
    let Some(offset) = offset else {
        return Ok(None);
    };
    let array = batch.column(offset);
    if array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<Float64Array>() {
        Ok(Some(array.value(row)))
    } else {
        bail!("numeric p-value columns must be Float64")
    }
}

fn optional_u32_value(
    batch: &RecordBatch,
    offset: Option<usize>,
    row: usize,
) -> Result<Option<u32>> {
    let Some(offset) = offset else {
        return Ok(None);
    };
    let array = batch.column(offset);
    if array.is_null(row) {
        return Ok(None);
    }
    if let Some(array) = array.as_any().downcast_ref::<UInt32Array>() {
        Ok(Some(array.value(row)))
    } else if let Some(array) = array.as_any().downcast_ref::<UInt64Array>() {
        Ok(Some(u32::try_from(array.value(row))?))
    } else if let Some(array) = array.as_any().downcast_ref::<Int64Array>() {
        Ok(Some(u32::try_from(array.value(row))?))
    } else {
        bail!("count columns must be UInt32, UInt64, or Int64")
    }
}
