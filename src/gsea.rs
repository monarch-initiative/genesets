//! Preranked gene set enrichment analysis (GSEA).
//!
//! Where Fisher exact and MGSA take a thresholded query set, GSEA takes every
//! measured gene with a score (a fold change, a signed -log10 p-value, a
//! t-statistic) and asks whether the members of each target set sit at the
//! top or the bottom of the ranking more than chance allows.
//!
//! The statistic follows Subramanian et al. (2005), "Gene set enrichment
//! analysis: a knowledge-based approach for interpreting genome-wide
//! expression profiles":
//!
//! - genes are sorted by score, highest first;
//! - a running sum walks down the list, stepping up by `|score|^weight / NR`
//!   at each set member (NR is the sum over members) and down by `1 / (N - k)`
//!   at every other gene;
//! - the enrichment score (ES) is the largest deviation of the running sum
//!   from zero, positive for sets near the top and negative near the bottom.
//!   `weight = 0` gives the unweighted Kolmogorov-Smirnov statistic and
//!   `weight = 1` is the GSEA default.
//!
//! Significance comes from gene-set permutation, as in GSEA preranked and
//! fgsea's simple mode: the null for a set of size `k` is the ES of random
//! `k`-gene subsets of the same ranking. Every target with the same in-ranking
//! size shares one null, so the cost scales with the number of distinct sizes
//! rather than the number of targets. Each null is drawn from its own random
//! stream, so results do not depend on the thread count.
//!
//! The p-value is computed within the null of the same sign as the observed
//! ES, the normalized enrichment score (NES) divides the ES by the mean of
//! that same-sign null, and `p_adjust` is Benjamini-Hochberg over the targets
//! tested against one ranking.

use std::collections::{BTreeMap, HashSet};

use anyhow::{Result, anyhow, bail};
use rayon::prelude::*;

use crate::{
    index::{GeneUniverse, NamedSet},
    rng::Rng,
};

/// One scored gene list, in input order.
#[derive(Clone, Debug, PartialEq)]
pub struct RankedList {
    pub id: String,
    pub genes: Vec<String>,
    pub scores: Vec<f64>,
}

impl RankedList {
    pub fn new(id: impl Into<String>, genes: Vec<String>, scores: Vec<f64>) -> Self {
        Self {
            id: id.into(),
            genes,
            scores,
        }
    }
}

/// A ranking resolved against a gene universe and sorted by score, highest
/// first. Ties keep input order.
#[derive(Clone, Debug)]
pub struct Ranking {
    pub id: String,
    /// Universe index of the gene at each rank.
    pub genes: Vec<usize>,
    /// Score at each rank, descending.
    pub scores: Vec<f64>,
    /// Rank of each universe gene, or `None` if the gene is not ranked.
    positions: Vec<Option<u32>>,
}

impl Ranking {
    /// Sorts `list` by score and indexes it. Every gene must already be in
    /// `universe`, appear once, and have a finite score.
    pub fn new(list: &RankedList, universe: &GeneUniverse) -> Result<Self> {
        if list.genes.len() != list.scores.len() {
            bail!(
                "ranking {:?} has {} genes but {} scores",
                list.id,
                list.genes.len(),
                list.scores.len()
            );
        }
        let mut seen = HashSet::new();
        let mut entries = Vec::with_capacity(list.genes.len());
        for (gene, score) in list.genes.iter().zip(&list.scores) {
            if !score.is_finite() {
                bail!("ranking {:?} has a non-finite score for {gene}", list.id);
            }
            if !seen.insert(gene.as_str()) {
                bail!("ranking {:?} lists {gene} more than once", list.id);
            }
            let index = universe
                .get(gene)
                .ok_or_else(|| anyhow!("ranked gene {gene} is not in the universe"))?;
            entries.push((index, *score));
        }
        entries.sort_by(|left, right| right.1.total_cmp(&left.1));

        let mut positions = vec![None; universe.len()];
        for (rank, (index, _)) in entries.iter().enumerate() {
            positions[*index] = Some(rank as u32);
        }
        Ok(Self {
            id: list.id.clone(),
            genes: entries.iter().map(|(index, _)| *index).collect(),
            scores: entries.iter().map(|(_, score)| *score).collect(),
            positions,
        })
    }

    pub fn len(&self) -> usize {
        self.genes.len()
    }

    pub fn is_empty(&self) -> bool {
        self.genes.is_empty()
    }

    /// Ranks of the members of `set` that are in this ranking, ascending.
    pub fn hit_positions(&self, set: &NamedSet) -> Vec<u32> {
        let mut hits: Vec<u32> = set
            .bits
            .ones()
            .into_iter()
            .filter_map(|index| self.positions.get(index).copied().flatten())
            .collect();
        hits.sort_unstable();
        hits
    }
}

#[derive(Clone, Debug)]
pub struct GseaOptions {
    /// Exponent on `|score|` when stepping up at a set member.
    pub weight: f64,
    /// Random gene sets drawn for each distinct set size.
    pub permutations: u32,
    /// Targets with fewer ranked members are not tested.
    pub min_size: u32,
    /// Targets with more ranked members are not tested.
    pub max_size: u32,
    pub seed: u64,
    pub include_leading_edge: bool,
    pub max_p_value: Option<f64>,
    pub max_p_adjust: Option<f64>,
}

impl Default for GseaOptions {
    fn default() -> Self {
        Self {
            weight: 1.0,
            permutations: 1000,
            min_size: 15,
            max_size: 500,
            seed: 123_456_789,
            include_leading_edge: false,
            max_p_value: None,
            max_p_adjust: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct GseaRow {
    pub query_index: usize,
    pub target_index: usize,
    /// Target members present in the ranking.
    pub set_size: u32,
    /// Target members in the universe, ranked or not.
    pub target_size: u32,
    pub ranked_size: u32,
    pub es: f64,
    pub nes: f64,
    pub p_value: f64,
    /// Benjamini-Hochberg over the targets tested against this ranking.
    pub p_adjust: f64,
    /// Members up to the ES peak (positive ES) or from the trough on
    /// (negative ES); these drive the score.
    pub leading_edge_size: u32,
    /// Rank, counted from 1, at which the running sum reaches the ES.
    pub peak_rank: u32,
    /// Universe indices of the leading edge, in rank order.
    pub leading_edge_indices: Option<Vec<usize>>,
}

/// Runs GSEA independently for each ranking. Rows are sorted by ranking, then
/// p-value, then by |NES| (largest first), then target id.
pub fn gsea_matrix(
    rankings: &[Ranking],
    targets: &[NamedSet],
    options: &GseaOptions,
) -> Vec<GseaRow> {
    rankings
        .par_iter()
        .enumerate()
        .flat_map(|(query_index, ranking)| gsea_ranking(query_index, ranking, targets, options))
        .collect()
}

pub fn gsea_one(ranking: &Ranking, targets: &[NamedSet], options: &GseaOptions) -> Vec<GseaRow> {
    gsea_ranking(0, ranking, targets, options)
}

struct Observed {
    target_index: usize,
    hits: Vec<u32>,
    score: RunningSumExtreme,
}

fn gsea_ranking(
    query_index: usize,
    ranking: &Ranking,
    targets: &[NamedSet],
    options: &GseaOptions,
) -> Vec<GseaRow> {
    let n = ranking.len();
    let weights: Vec<f64> = ranking
        .scores
        .iter()
        .map(|score| step_weight(*score, options.weight))
        .collect();

    let observed: Vec<Observed> = targets
        .par_iter()
        .enumerate()
        .filter_map(|(target_index, target)| {
            let hits = ranking.hit_positions(target);
            let k = hits.len();
            if k == 0 || k >= n || k < options.min_size as usize || k > options.max_size as usize {
                return None;
            }
            let score = running_sum_extreme(&hits, &weights, n);
            Some(Observed {
                target_index,
                hits,
                score,
            })
        })
        .collect();

    let mut sizes: Vec<usize> = observed.iter().map(|item| item.hits.len()).collect();
    sizes.sort_unstable();
    sizes.dedup();
    let query_seed = options.seed ^ (query_index as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
    let nulls: BTreeMap<usize, NullDistribution> = sizes
        .par_iter()
        .map(|k| {
            let seed = query_seed ^ (*k as u64).wrapping_mul(0xD1B5_4A32_D192_ED03);
            (
                *k,
                NullDistribution::sample(*k, &weights, options.permutations, seed),
            )
        })
        .collect();

    let mut rows: Vec<GseaRow> = observed
        .into_iter()
        .map(|item| {
            let null = &nulls[&item.hits.len()];
            let es = item.score.es;
            let (p_value, nes) = null.test(es);
            let (leading_edge, peak_rank) = if es >= 0.0 {
                let end = item.score.max_hit + 1;
                (&item.hits[..end], item.hits[item.score.max_hit] + 1)
            } else {
                // The trough sits on the miss just before this hit.
                let start = item.score.min_hit;
                (&item.hits[start..], item.hits[start])
            };
            GseaRow {
                query_index,
                target_index: item.target_index,
                set_size: item.hits.len() as u32,
                target_size: targets[item.target_index].bits.count_ones(),
                ranked_size: n as u32,
                es,
                nes,
                p_value,
                p_adjust: f64::NAN,
                leading_edge_size: leading_edge.len() as u32,
                peak_rank,
                leading_edge_indices: options.include_leading_edge.then(|| {
                    leading_edge
                        .iter()
                        .map(|rank| ranking.genes[*rank as usize])
                        .collect()
                }),
            }
        })
        .collect();

    benjamini_hochberg(&mut rows);
    rows.retain(|row| {
        options
            .max_p_value
            .is_none_or(|cutoff| row.p_value <= cutoff)
            && options
                .max_p_adjust
                .is_none_or(|cutoff| row.p_adjust <= cutoff)
    });
    rows.sort_by(|left, right| {
        left.p_value
            .total_cmp(&right.p_value)
            .then_with(|| right.nes.abs().total_cmp(&left.nes.abs()))
            .then_with(|| {
                targets[left.target_index]
                    .id
                    .cmp(&targets[right.target_index].id)
            })
    });
    rows
}

fn step_weight(score: f64, weight: f64) -> f64 {
    if weight == 0.0 {
        1.0
    } else {
        score.abs().powf(weight)
    }
}

#[derive(Copy, Clone, Debug, PartialEq)]
struct RunningSumExtreme {
    es: f64,
    /// Hit at which the running sum peaks.
    max_hit: usize,
    /// Hit just before which the running sum bottoms out.
    min_hit: usize,
}

/// Walks the running sum over the sorted `hits`. The sum only rises at a hit,
/// so its maximum is right after some hit and its minimum right before one
/// (or at the end, where it returns to zero).
fn running_sum_extreme(hits: &[u32], weights: &[f64], n: usize) -> RunningSumExtreme {
    let k = hits.len();
    let hit_total: f64 = hits.iter().map(|rank| weights[*rank as usize]).sum();
    // With all member scores at zero there is nothing to weight by, so every
    // member steps up equally, as in the unweighted statistic.
    let uniform = hit_total <= 0.0;
    let miss_step = 1.0 / (n - k) as f64;

    let mut hit_sum = 0.0;
    let mut max_dev = 0.0;
    let mut min_dev = 0.0;
    let mut max_hit = 0;
    let mut min_hit = 0;
    for (j, rank) in hits.iter().enumerate() {
        let misses = (*rank as usize - j) as f64 * miss_step;
        let before = hit_sum - misses;
        if before < min_dev {
            min_dev = before;
            min_hit = j;
        }
        hit_sum += if uniform {
            1.0 / k as f64
        } else {
            weights[*rank as usize] / hit_total
        };
        let after = hit_sum - misses;
        if after > max_dev {
            max_dev = after;
            max_hit = j;
        }
    }
    RunningSumExtreme {
        es: if max_dev >= -min_dev {
            max_dev
        } else {
            min_dev
        },
        max_hit,
        min_hit,
    }
}

/// Permutation null for one set size, split by sign.
struct NullDistribution {
    /// Non-negative null ES values, ascending.
    positive: Vec<f64>,
    /// Negative null ES values, ascending.
    negative: Vec<f64>,
}

impl NullDistribution {
    fn sample(k: usize, weights: &[f64], permutations: u32, seed: u64) -> Self {
        let n = weights.len();
        let mut rng = Rng::new(seed);
        let mut deck: Vec<u32> = (0..n as u32).collect();
        let mut swaps = Vec::with_capacity(k);
        let mut hits = Vec::with_capacity(k);
        let mut positive = Vec::new();
        let mut negative = Vec::new();
        for _ in 0..permutations {
            // Partial Fisher-Yates, undone afterwards so each draw is O(k).
            swaps.clear();
            for i in 0..k {
                let j = i + rng.index(n - i);
                deck.swap(i, j);
                swaps.push(j);
            }
            hits.clear();
            hits.extend_from_slice(&deck[..k]);
            for (i, j) in swaps.iter().enumerate().rev() {
                deck.swap(i, *j);
            }
            hits.sort_unstable();
            let es = running_sum_extreme(&hits, weights, n).es;
            if es >= 0.0 {
                positive.push(es);
            } else {
                negative.push(es);
            }
        }
        positive.sort_by(f64::total_cmp);
        negative.sort_by(f64::total_cmp);
        Self { positive, negative }
    }

    /// Returns the permutation p-value and the NES for an observed ES.
    fn test(&self, es: f64) -> (f64, f64) {
        if es >= 0.0 {
            let at_least = self.positive.len() - self.positive.partition_point(|null| *null < es);
            let p_value = (at_least + 1) as f64 / (self.positive.len() + 1) as f64;
            (p_value, es / mean(&self.positive))
        } else {
            let at_most = self.negative.partition_point(|null| *null <= es);
            let p_value = (at_most + 1) as f64 / (self.negative.len() + 1) as f64;
            (p_value, es / mean(&self.negative).abs())
        }
    }
}

fn mean(values: &[f64]) -> f64 {
    if values.is_empty() {
        f64::NAN
    } else {
        values.iter().sum::<f64>() / values.len() as f64
    }
}

fn benjamini_hochberg(rows: &mut [GseaRow]) {
    let m = rows.len();
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|left, right| rows[*right].p_value.total_cmp(&rows[*left].p_value));
    let mut running_min: f64 = 1.0;
    for (offset, index) in order.into_iter().enumerate() {
        let rank = m - offset;
        let adjusted = rows[index].p_value * m as f64 / rank as f64;
        running_min = running_min.min(adjusted);
        rows[index].p_adjust = running_min;
    }
}

#[cfg(test)]
mod tests {
    use super::{
        GseaOptions, RankedList, Ranking, benjamini_hochberg, gsea_one, running_sum_extreme,
    };
    use crate::index::{GeneUniverse, RawSet, build_flat_sets};

    fn approx(left: f64, right: f64) -> bool {
        (left - right).abs() < 1e-12
    }

    #[test]
    fn running_sum_matches_hand_computation() {
        // Ten genes, scores 10..1, members at ranks 0, 2 and 7.
        let weights: Vec<f64> = (1..=10).rev().map(f64::from).collect();
        let hits = [0, 2, 7];
        let extreme = running_sum_extreme(&hits, &weights, weights.len());
        // NR = 10 + 8 + 3; peak after rank 2: (10 + 8) / 21 - 1 / 7.
        assert!(approx(extreme.es, 18.0 / 21.0 - 1.0 / 7.0));
        assert_eq!(extreme.max_hit, 1);

        // Unweighted, members at the bottom: the trough is right before them.
        let flat = vec![1.0; 10];
        let extreme = running_sum_extreme(&[7, 8, 9], &flat, 10);
        assert!(approx(extreme.es, -1.0));
        assert_eq!(extreme.min_hit, 0);
    }

    fn fixture() -> (GeneUniverse, Ranking, Vec<crate::index::NamedSet>) {
        let genes: Vec<String> = (0..200).map(|i| format!("g{i:03}")).collect();
        let scores: Vec<f64> = (0..200).map(|i| 100.0 - i as f64).collect();
        let mut universe = GeneUniverse::new();
        universe.add_genes(&genes);
        let list = RankedList::new("contrast", genes.clone(), scores);
        let ranking = Ranking::new(&list, &universe).unwrap();
        let pick = |range: std::ops::Range<usize>, step: usize| -> Vec<String> {
            range.step_by(step).map(|i| genes[i].clone()).collect()
        };
        let raw = vec![
            RawSet::new("top", None, pick(0..40, 2)),
            RawSet::new("bottom", None, pick(160..200, 2)),
            RawSet::new("spread", None, pick(0..200, 10)),
        ];
        let sets = build_flat_sets(&raw, &universe);
        (universe, ranking, sets)
    }

    #[test]
    fn separates_top_bottom_and_spread_sets() {
        let (_, ranking, sets) = fixture();
        let options = GseaOptions {
            permutations: 2000,
            include_leading_edge: true,
            ..GseaOptions::default()
        };
        let rows = gsea_one(&ranking, &sets, &options);
        assert_eq!(rows.len(), 3);
        let row = |id: &str| {
            rows.iter()
                .find(|row| sets[row.target_index].id == id)
                .unwrap()
        };

        let top = row("top");
        assert!(top.es > 0.8 && top.nes > 1.0);
        assert!(top.p_value < 0.005 && top.p_adjust < 0.01);
        assert_eq!(top.leading_edge_size, 20);
        assert_eq!(top.peak_rank, 39);

        let bottom = row("bottom");
        assert!(bottom.es < -0.8 && bottom.nes < -1.0);
        assert!(bottom.p_value < 0.005);
        assert_eq!(bottom.leading_edge_size, 20);
        assert_eq!(bottom.peak_rank, 160);

        assert!(row("spread").p_value > 0.1);
        assert_eq!(
            rows.last().map(|row| sets[row.target_index].id.as_str()),
            Some("spread")
        );
    }

    #[test]
    fn is_reproducible_and_respects_size_limits() {
        let (_, ranking, sets) = fixture();
        let options = GseaOptions {
            permutations: 200,
            ..GseaOptions::default()
        };
        let first = gsea_one(&ranking, &sets, &options);
        let second = gsea_one(&ranking, &sets, &options);
        let summary =
            |rows: &[super::GseaRow]| rows.iter().map(|row| row.p_value).collect::<Vec<_>>();
        assert_eq!(summary(&first), summary(&second));

        let limited = GseaOptions {
            min_size: 21,
            ..options
        };
        assert!(gsea_one(&ranking, &sets, &limited).is_empty());
    }

    #[test]
    fn rejects_duplicate_and_non_finite_scores() {
        let mut universe = GeneUniverse::new();
        let genes = vec!["a".to_owned(), "a".to_owned()];
        universe.add_genes(&genes);
        let duplicate = RankedList::new("r", genes, vec![1.0, 2.0]);
        assert!(Ranking::new(&duplicate, &universe).is_err());
        let nan = RankedList::new("r", vec!["a".to_owned()], vec![f64::NAN]);
        assert!(Ranking::new(&nan, &universe).is_err());
    }

    #[test]
    fn adjusts_with_benjamini_hochberg() {
        let (_, ranking, sets) = fixture();
        let mut rows = gsea_one(&ranking, &sets, &GseaOptions::default());
        for (row, p_value) in rows.iter_mut().zip([0.01, 0.04, 0.03]) {
            row.p_value = p_value;
        }
        benjamini_hochberg(&mut rows);
        let adjusted: Vec<f64> = rows.iter().map(|row| row.p_adjust).collect();
        assert!(approx(adjusted[0], 0.03));
        assert!(approx(adjusted[1], 0.04));
        assert!(approx(adjusted[2], 0.04));
    }
}
