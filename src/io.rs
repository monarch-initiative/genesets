use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader},
    path::{Path, PathBuf},
};

use anyhow::{Context, Result, bail};
use clap::ValueEnum;
use csv::{ReaderBuilder, Trim};
use serde::Deserialize;

use crate::index::RawSet;

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[value(rename_all = "kebab-case")]
pub enum SetFormat {
    Auto,
    List,
    Pairwise,
    GeneTerm,
    Gmt,
    Gmx,
    GmxDesc,
}

pub fn read_name_pairs(path: &Path) -> Result<Vec<(String, String)>> {
    let records = read_records(path)?;
    let mut pairs = Vec::new();
    for (row_index, record) in records.into_iter().enumerate() {
        if record.len() < 2 {
            continue;
        }
        if row_index == 0
            && looks_like_pair_header(
                &record,
                &[
                    "term",
                    "termid",
                    "curie",
                    "id",
                    "gene",
                    "geneid",
                    "genecurie",
                ],
                &["name", "label", "termname", "symbol"],
            )
        {
            continue;
        }
        if !record[0].is_empty() {
            pairs.push((record[0].clone(), record[1].clone()));
        }
    }
    Ok(pairs)
}

pub fn read_name_map(path: &Path) -> Result<HashMap<String, String>> {
    Ok(read_name_pairs(path)?.into_iter().collect())
}

pub fn read_annotation_pairs(path: &Path) -> Result<Vec<(String, String)>> {
    read_pair_table(
        path,
        &["gene", "geneid", "genecurie", "curie"],
        &["term", "termid", "termcurie", "set", "setid"],
    )
}

pub fn read_closure_pairs(path: &Path) -> Result<Vec<(String, String)>> {
    read_pair_table(
        path,
        &["child", "subject", "descendant", "term", "termid"],
        &["ancestor", "parent", "object", "superterm", "superclass"],
    )
}

pub fn read_gene_list(path: &Path) -> Result<Vec<String>> {
    let records = read_records(path)?;
    let mut genes = Vec::new();
    for (row_index, record) in records.into_iter().enumerate() {
        if record.is_empty() || record[0].is_empty() {
            continue;
        }
        if row_index == 0 && is_one_column_header(&record[0]) {
            continue;
        }
        genes.push(record[0].clone());
    }
    Ok(genes)
}

pub fn read_sets(path: &Path, format: SetFormat, default_name: &str) -> Result<Vec<RawSet>> {
    let format = if format == SetFormat::Auto {
        detect_set_format(path)?
    } else {
        format
    };

    match format {
        SetFormat::Auto => unreachable!("auto set format is resolved before parsing"),
        SetFormat::List => read_list_set(path, default_name),
        SetFormat::Pairwise => read_pairwise_sets(path, PairwiseOrientation::SetGene),
        SetFormat::GeneTerm => read_pairwise_sets(path, PairwiseOrientation::GeneSet),
        SetFormat::Gmt => read_gmt(path),
        SetFormat::Gmx => read_gmx(path, false),
        SetFormat::GmxDesc => read_gmx(path, true),
    }
}

fn read_pair_table(
    path: &Path,
    first_header_words: &[&str],
    second_header_words: &[&str],
) -> Result<Vec<(String, String)>> {
    let records = read_records(path)?;
    let mut pairs = Vec::new();
    for (row_index, record) in records.into_iter().enumerate() {
        if record.len() < 2 {
            continue;
        }
        if row_index == 0
            && looks_like_pair_header(&record, first_header_words, second_header_words)
        {
            continue;
        }
        if !record[0].is_empty() && !record[1].is_empty() {
            pairs.push((record[0].clone(), record[1].clone()));
        }
    }
    Ok(pairs)
}

fn read_list_set(path: &Path, default_name: &str) -> Result<Vec<RawSet>> {
    Ok(vec![RawSet::new(default_name, None, read_gene_list(path)?)])
}

#[derive(Copy, Clone, Debug, Eq, PartialEq)]
enum PairwiseOrientation {
    SetGene,
    GeneSet,
}

fn read_pairwise_sets(path: &Path, orientation: PairwiseOrientation) -> Result<Vec<RawSet>> {
    let records = read_records(path)?;
    let mut sets = Vec::<RawSet>::new();
    let mut offsets = HashMap::<String, usize>::new();

    for (row_index, record) in records.into_iter().enumerate() {
        if record.len() < 2 {
            continue;
        }
        if row_index == 0 {
            let is_header = match orientation {
                PairwiseOrientation::SetGene => looks_like_pair_header(
                    &record,
                    &["set", "setid", "term", "termid"],
                    &["gene", "geneid", "genecurie", "curie"],
                ),
                PairwiseOrientation::GeneSet => looks_like_pair_header(
                    &record,
                    &["gene", "geneid", "genecurie", "curie"],
                    &["set", "setid", "term", "termid"],
                ),
            };
            if is_header {
                continue;
            }
        }

        let (set_id, gene) = match orientation {
            PairwiseOrientation::SetGene => (&record[0], &record[1]),
            PairwiseOrientation::GeneSet => (&record[1], &record[0]),
        };
        if set_id.is_empty() || gene.is_empty() {
            continue;
        }

        let set_index = if let Some(index) = offsets.get(set_id) {
            *index
        } else {
            let index = sets.len();
            offsets.insert(set_id.clone(), index);
            sets.push(RawSet::new(set_id.clone(), None, Vec::new()));
            index
        };
        sets[set_index].genes.push(gene.clone());
    }

    Ok(sets)
}

fn read_gmt(path: &Path) -> Result<Vec<RawSet>> {
    let records = read_records(path)?;
    let mut sets = Vec::new();
    for record in records {
        if record.len() < 2 || record[0].is_empty() {
            continue;
        }
        let has_description = record.len() >= 3;
        let name = has_description
            .then(|| record[1].clone())
            .filter(|name| !name.is_empty());
        let genes_start = if has_description { 2 } else { 1 };
        let genes = record[genes_start..]
            .iter()
            .filter(|gene| !gene.is_empty())
            .cloned()
            .collect();
        sets.push(RawSet::new(record[0].clone(), name, genes));
    }
    Ok(sets)
}

fn read_gmx(path: &Path, has_description_row: bool) -> Result<Vec<RawSet>> {
    let records = read_records(path)?;
    if records.is_empty() {
        bail!("GMX file {} is empty", path.display());
    }

    let set_ids = &records[0];
    let mut sets: Vec<RawSet> = set_ids
        .iter()
        .filter(|id| !id.is_empty())
        .map(|id| RawSet::new(id.clone(), None, Vec::new()))
        .collect();
    if sets.is_empty() {
        bail!(
            "GMX file {} does not contain any set ids in the first row",
            path.display()
        );
    }

    let data_start = if has_description_row {
        if let Some(description_row) = records.get(1) {
            for (set_index, description) in description_row.iter().enumerate().take(sets.len()) {
                if !description.is_empty() {
                    sets[set_index].name = Some(description.clone());
                }
            }
        }
        2
    } else {
        1
    };

    for record in records.iter().skip(data_start) {
        for (set_index, gene) in record.iter().enumerate().take(sets.len()) {
            if !gene.is_empty() {
                sets[set_index].genes.push(gene.clone());
            }
        }
    }

    Ok(sets)
}

fn detect_set_format(path: &Path) -> Result<SetFormat> {
    match path.extension().and_then(|extension| extension.to_str()) {
        Some(extension) if extension.eq_ignore_ascii_case("gmt") => return Ok(SetFormat::Gmt),
        Some(extension) if extension.eq_ignore_ascii_case("gmx") => return Ok(SetFormat::Gmx),
        _ => {}
    }

    let records = read_records(path)?;
    let Some(first) = records.first() else {
        bail!("cannot detect format for empty set file {}", path.display());
    };
    Ok(match first.len() {
        0 | 1 => SetFormat::List,
        2 => SetFormat::Pairwise,
        _ => SetFormat::Gmt,
    })
}

fn read_records(path: &Path) -> Result<Vec<Vec<String>>> {
    let delimiter = detect_delimiter(path)?;
    let mut reader = ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(false)
        .flexible(true)
        .comment(Some(b'#'))
        .trim(Trim::All)
        .from_path(path)
        .with_context(|| format!("failed to open {}", path.display()))?;

    let mut records = Vec::new();
    for record in reader.records() {
        let record = record.with_context(|| format!("failed to parse {}", path.display()))?;
        let fields: Vec<String> = record.iter().map(|field| field.trim().to_owned()).collect();
        if fields.iter().all(|field| field.is_empty()) {
            continue;
        }
        records.push(fields);
    }
    Ok(records)
}

fn detect_delimiter(path: &Path) -> Result<u8> {
    if path
        .extension()
        .and_then(|extension| extension.to_str())
        .is_some_and(|extension| extension.eq_ignore_ascii_case("csv"))
    {
        return Ok(b',');
    }

    let file = File::open(path).with_context(|| format!("failed to open {}", path.display()))?;
    for line in BufReader::new(file).lines() {
        let line = line.with_context(|| format!("failed to read {}", path.display()))?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        if line.contains(',') && !line.contains('\t') {
            return Ok(b',');
        }
        return Ok(b'\t');
    }
    Ok(b'\t')
}

fn looks_like_pair_header(record: &[String], first_words: &[&str], second_words: &[&str]) -> bool {
    if record.len() < 2 {
        return false;
    }
    let first = normalize_header(&record[0]);
    let second = normalize_header(&record[1]);
    first_words.contains(&first.as_str()) && second_words.contains(&second.as_str())
}

fn is_one_column_header(value: &str) -> bool {
    matches!(
        normalize_header(value).as_str(),
        "gene" | "geneid" | "genecurie" | "curie" | "id"
    )
}

fn normalize_header(value: &str) -> String {
    value
        .trim()
        .to_ascii_lowercase()
        .chars()
        .filter(|character| !matches!(character, '_' | '-' | ' '))
        .collect()
}

pub fn default_set_name(path: &Path) -> String {
    path.file_stem()
        .and_then(|stem| stem.to_str())
        .map(str::to_owned)
        .unwrap_or_else(|| PathBuf::from(path).display().to_string())
}

#[cfg(test)]
mod tests {
    use std::{fs, path::PathBuf};

    use super::{SetFormat, read_sets};

    fn temp_file(name: &str, contents: &str) -> PathBuf {
        let path =
            std::env::temp_dir().join(format!("genesets-rs-{}-{}", std::process::id(), name));
        fs::write(&path, contents).unwrap();
        path
    }

    #[test]
    fn reads_pairwise_sets() {
        let path = temp_file("pairwise.tsv", "set\tgene\ns1\tg1\ns1\tg2\ns2\tg2\n");
        let sets = read_sets(&path, SetFormat::Pairwise, "sample").unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(sets.len(), 2);
        assert_eq!(sets[0].id, "s1");
        assert_eq!(sets[0].genes, vec!["g1", "g2"]);
    }

    #[test]
    fn reads_simple_gmx_without_description_row() {
        let path = temp_file("sets.gmx", "s1\ts2\ng1\tg3\ng2\t\n");
        let sets = read_sets(&path, SetFormat::Gmx, "sample").unwrap();
        fs::remove_file(path).unwrap();

        assert_eq!(sets.len(), 2);
        assert_eq!(sets[0].genes, vec!["g1", "g2"]);
        assert_eq!(sets[1].genes, vec!["g3"]);
    }
}
