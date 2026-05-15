use clap::ValueEnum;
use rayon::prelude::*;
use serde::Deserialize;

use crate::{bitset::DenseBitSet, fisher::FisherCache, index::NamedSet};

#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[value(rename_all = "kebab-case")]
pub enum Correction {
    Bonferroni,
    None,
}

#[derive(Clone, Debug)]
pub struct EnrichmentOptions {
    pub min_overlap: u32,
    pub correction: Correction,
    pub include_overlap_genes: bool,
    pub max_p_value: Option<f64>,
    pub max_p_adjust: Option<f64>,
}

impl Default for EnrichmentOptions {
    fn default() -> Self {
        Self {
            min_overlap: 1,
            correction: Correction::Bonferroni,
            include_overlap_genes: false,
            max_p_value: None,
            max_p_adjust: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct EnrichmentRow {
    pub query_index: usize,
    pub target_index: usize,
    pub overlap: u32,
    pub query_size: u32,
    pub target_size: u32,
    pub background_size: u32,
    pub p_value: f64,
    pub p_adjust: f64,
    pub overlap_indices: Option<Vec<usize>>,
}

pub fn enrich_one(
    query: &NamedSet,
    targets: &[NamedSet],
    background: &DenseBitSet,
    options: &EnrichmentOptions,
) -> Vec<EnrichmentRow> {
    let rows = enrich_matrix(std::slice::from_ref(query), targets, background, options);
    rows.into_iter()
        .map(|mut row| {
            row.query_index = 0;
            row
        })
        .collect()
}

pub fn enrich_matrix(
    queries: &[NamedSet],
    targets: &[NamedSet],
    background: &DenseBitSet,
    options: &EnrichmentOptions,
) -> Vec<EnrichmentRow> {
    let background_size = background.count_ones();
    let query_sizes: Vec<u32> = queries
        .iter()
        .map(|query| query.bits.and_count(background))
        .collect();
    let target_sizes: Vec<u32> = targets
        .iter()
        .map(|target| target.bits.and_count(background))
        .collect();

    let non_empty_queries = query_sizes.iter().filter(|size| **size > 0).count();
    let non_empty_targets = target_sizes.iter().filter(|size| **size > 0).count();
    let test_count = (non_empty_queries * non_empty_targets).max(1) as f64;
    let fisher = FisherCache::new(background_size as usize);

    let mut rows: Vec<EnrichmentRow> = queries
        .par_iter()
        .enumerate()
        .flat_map(|(query_index, query)| {
            targets
                .iter()
                .enumerate()
                .filter_map(|(target_index, target)| {
                    compute_row(
                        query_index,
                        query,
                        query_sizes[query_index],
                        target_index,
                        target,
                        target_sizes[target_index],
                        background,
                        background_size,
                        test_count,
                        &fisher,
                        options,
                    )
                })
                .collect::<Vec<_>>()
        })
        .collect();

    if queries.len() == 1 {
        rows.sort_by(|left, right| {
            left.p_value
                .total_cmp(&right.p_value)
                .then_with(|| right.overlap.cmp(&left.overlap))
                .then_with(|| {
                    targets[left.target_index]
                        .id
                        .cmp(&targets[right.target_index].id)
                })
        });
    } else {
        rows.sort_by(|left, right| {
            left.query_index
                .cmp(&right.query_index)
                .then_with(|| left.p_value.total_cmp(&right.p_value))
                .then_with(|| left.target_index.cmp(&right.target_index))
        });
    }

    rows
}

#[allow(clippy::too_many_arguments)]
fn compute_row(
    query_index: usize,
    query: &NamedSet,
    query_size: u32,
    target_index: usize,
    target: &NamedSet,
    target_size: u32,
    background: &DenseBitSet,
    background_size: u32,
    test_count: f64,
    fisher: &FisherCache,
    options: &EnrichmentOptions,
) -> Option<EnrichmentRow> {
    let overlap = query.bits.and3_count(&target.bits, background);
    if overlap < options.min_overlap {
        return None;
    }

    let p_value = if background_size == 0 || query_size == 0 || target_size == 0 {
        1.0
    } else {
        fisher.right_tail(background_size, target_size, query_size, overlap)
    };
    let p_adjust = match options.correction {
        Correction::Bonferroni => (p_value * test_count).min(1.0),
        Correction::None => p_value,
    };
    if options
        .max_p_value
        .is_some_and(|max_p_value| p_value > max_p_value)
        || options
            .max_p_adjust
            .is_some_and(|max_p_adjust| p_adjust > max_p_adjust)
    {
        return None;
    }
    let overlap_indices = options
        .include_overlap_genes
        .then(|| query.bits.intersection_indices3(&target.bits, background));

    Some(EnrichmentRow {
        query_index,
        target_index,
        overlap,
        query_size,
        target_size,
        background_size,
        p_value,
        p_adjust,
        overlap_indices,
    })
}

#[cfg(test)]
mod tests {
    use crate::{
        enrichment::{Correction, EnrichmentOptions, enrich_one},
        index::{GeneUniverse, RawSet, build_flat_sets, build_gene_list_bitset},
    };

    #[test]
    fn enriches_single_query_against_flat_targets() {
        let mut universe = GeneUniverse::new();
        for gene in ["g1", "g2", "g3", "g4", "g5"] {
            universe.add_gene(gene);
        }

        let query = build_flat_sets(
            &[RawSet::new(
                "sample",
                None,
                vec!["g1".into(), "g2".into(), "g5".into()],
            )],
            &universe,
        )
        .remove(0);
        let targets = build_flat_sets(
            &[
                RawSet::new("child", None, vec!["g1".into(), "g2".into()]),
                RawSet::new("parent", None, vec!["g1".into(), "g2".into(), "g3".into()]),
            ],
            &universe,
        );
        let background = build_gene_list_bitset(
            &[
                "g1".into(),
                "g2".into(),
                "g3".into(),
                "g4".into(),
                "g5".into(),
            ],
            &universe,
        );
        let rows = enrich_one(
            &query,
            &targets,
            &background,
            &EnrichmentOptions {
                min_overlap: 1,
                correction: Correction::Bonferroni,
                include_overlap_genes: true,
                max_p_value: None,
                max_p_adjust: None,
            },
        );

        assert_eq!(rows.len(), 2);
        assert_eq!(targets[rows[0].target_index].id, "child");
        assert!((rows[0].p_value - 0.3).abs() < 1e-12);
        assert!((rows[0].p_adjust - 0.6).abs() < 1e-12);
        assert_eq!(rows[0].overlap_indices.as_ref().unwrap().len(), 2);
    }

    #[test]
    fn filters_by_adjusted_p_value() {
        let mut universe = GeneUniverse::new();
        for gene in ["g1", "g2", "g3", "g4", "g5"] {
            universe.add_gene(gene);
        }

        let query = build_flat_sets(
            &[RawSet::new(
                "sample",
                None,
                vec!["g1".into(), "g2".into(), "g5".into()],
            )],
            &universe,
        )
        .remove(0);
        let targets = build_flat_sets(
            &[
                RawSet::new("child", None, vec!["g1".into(), "g2".into()]),
                RawSet::new("parent", None, vec!["g1".into(), "g2".into(), "g3".into()]),
            ],
            &universe,
        );
        let background = build_gene_list_bitset(
            &[
                "g1".into(),
                "g2".into(),
                "g3".into(),
                "g4".into(),
                "g5".into(),
            ],
            &universe,
        );
        let rows = enrich_one(
            &query,
            &targets,
            &background,
            &EnrichmentOptions {
                min_overlap: 1,
                correction: Correction::Bonferroni,
                include_overlap_genes: false,
                max_p_value: None,
                max_p_adjust: Some(0.75),
            },
        );

        assert_eq!(rows.len(), 1);
        assert_eq!(targets[rows[0].target_index].id, "child");
    }
}
