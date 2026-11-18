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

use std::{
    collections::{BTreeMap, HashSet},
    sync::Arc,
};

use anyhow::{Result, anyhow, bail};
use clap::ValueEnum;
use rayon::prelude::*;
use serde::Deserialize;

use crate::{
    index::{GeneUniverse, NamedSet, TermGraph},
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

/// How a target set is scored against a ranking.
#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[value(rename_all = "kebab-case")]
pub enum Statistic {
    /// Weighted Kolmogorov-Smirnov running sum (GSEA) with permutation
    /// p-values.
    Ks,
    /// Wilcoxon rank-sum (Mann-Whitney U) of members against non-members,
    /// with a tie-corrected normal approximation.
    Wilcoxon,
}

/// How the ontology graph is used while scoring.
#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum, Deserialize)]
#[serde(rename_all = "kebab-case")]
#[value(rename_all = "kebab-case")]
pub enum Algorithm {
    /// Score every term on its own.
    Classic,
    /// Remove the genes of descendants with p at or below the elim cutoff.
    Elim,
    /// Remove a child's genes from its parent when the child is at least as
    /// significant as the parent.
    Weight01,
}

#[derive(Clone, Debug)]
pub struct GseaOptions {
    pub statistic: Statistic,
    pub algorithm: Algorithm,
    /// Raw p-value at or below which `elim` removes a term's genes from its
    /// ancestors.
    pub elim_cutoff: f64,
    /// Exponent on `|score|` when stepping up at a set member.
    pub weight: f64,
    /// Random gene sets drawn for each KS null.
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
            statistic: Statistic::Ks,
            algorithm: Algorithm::Classic,
            elim_cutoff: 0.01,
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
    /// Target members scored: those in the ranking, less any eliminated.
    pub set_size: u32,
    /// Ranked members taken out by `elim` or `weight01`, along with the rest
    /// of the ranking's copies of them.
    pub eliminated: u32,
    /// Target members in the universe, ranked or not.
    pub target_size: u32,
    pub ranked_size: u32,
    /// KS: the enrichment score. Wilcoxon: the rank-biserial correlation.
    /// Both run from -1 to 1 and are positive toward high scores.
    pub es: f64,
    /// KS: ES over the mean same-sign null ES. Wilcoxon: the z-score.
    pub nes: f64,
    pub p_value: f64,
    /// Benjamini-Hochberg over the targets tested against this ranking.
    pub p_adjust: f64,
    /// Members up to the running-sum peak (positive) or from its trough on
    /// (negative). Wilcoxon uses the unweighted running sum.
    pub leading_edge_size: u32,
    /// Rank, counted from 1, of the leading-edge gene nearest the middle of
    /// the ranking: its last gene when positive, its first when negative.
    pub peak_rank: u32,
    /// Universe indices of the leading edge, in rank order.
    pub leading_edge_indices: Option<Vec<usize>>,
}

/// Runs GSEA independently for each ranking. Rows are sorted by ranking, then
/// p-value, then by |NES| (largest first), then target id.
///
/// `elim` and `weight01` need the term graph; `classic` ignores it.
pub fn gsea_matrix(
    rankings: &[Ranking],
    targets: &[NamedSet],
    graph: Option<&TermGraph>,
    options: &GseaOptions,
) -> Result<Vec<GseaRow>> {
    if options.algorithm != Algorithm::Classic && graph.is_none() {
        bail!("the elim and weight01 algorithms need an ontology closure");
    }
    Ok(rankings
        .par_iter()
        .enumerate()
        .flat_map(|(query_index, ranking)| {
            gsea_ranking(query_index, ranking, targets, graph, options)
        })
        .collect())
}

pub fn gsea_one(
    ranking: &Ranking,
    targets: &[NamedSet],
    graph: Option<&TermGraph>,
    options: &GseaOptions,
) -> Result<Vec<GseaRow>> {
    gsea_matrix(std::slice::from_ref(ranking), targets, graph, options)
}

fn gsea_ranking(
    query_index: usize,
    ranking: &Ranking,
    targets: &[NamedSet],
    graph: Option<&TermGraph>,
    options: &GseaOptions,
) -> Vec<GseaRow> {
    let n = ranking.len();
    let hits: Vec<Option<Vec<u32>>> = targets
        .par_iter()
        .map(|target| {
            let hits = ranking.hit_positions(target);
            let k = hits.len();
            let tested =
                k > 0 && k < n && k >= options.min_size as usize && k <= options.max_size as usize;
            tested.then_some(hits)
        })
        .collect();

    let scorer = Scorer::new(query_index, ranking, &hits, options);
    let scores: Vec<Option<TermScore>> = match (options.algorithm, graph) {
        (Algorithm::Classic, _) | (_, None) => hits
            .par_iter()
            .enumerate()
            .map(|(target_index, hits)| {
                hits.as_ref()
                    .map(|hits| scorer.score(target_index, 0, hits, &[]))
            })
            .collect(),
        (algorithm, Some(graph)) => score_with_graph(&scorer, &hits, graph, algorithm, options),
    };

    let mut nulls = Vec::new();
    let mut rows: Vec<GseaRow> = scores
        .into_iter()
        .enumerate()
        .filter_map(|(target_index, score)| {
            let score = score?;
            nulls.push(score.null.clone());
            let ranked = hits[target_index].as_ref().map_or(0, Vec::len) as u32;
            Some(GseaRow {
                query_index,
                target_index,
                set_size: score.set_size,
                eliminated: ranked - score.set_size,
                target_size: targets[target_index].bits.count_ones(),
                ranked_size: n as u32,
                es: score.es,
                nes: score.nes,
                p_value: score.p_value,
                p_adjust: f64::NAN,
                leading_edge_size: score.leading_edge.len() as u32,
                peak_rank: score.peak_rank,
                leading_edge_indices: options.include_leading_edge.then(|| {
                    score
                        .leading_edge
                        .iter()
                        .map(|rank| ranking.genes[*rank as usize])
                        .collect()
                }),
            })
        })
        .collect();

    match options.statistic {
        Statistic::Ks => gsea_fdr(&mut rows, &nulls, &scorer.nulls),
        Statistic::Wilcoxon => benjamini_hochberg(&mut rows),
    }
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

/// A term's graph-aware result, kept for its ancestors.
struct GraphTerm {
    score: TermScore,
    /// Ranks this term asks its ancestors to remove.
    passed_up: Vec<u32>,
}

/// Scores tested terms deepest first, so every child is final before its
/// parent. A child has strictly more ancestors than its parent, so terms with
/// equal ancestor counts are never related and are scored in parallel.
///
/// `elim` follows Alexa et al. (2006) as topGO implements it: a term is scored
/// without the genes of descendants whose p-value is at or below the cutoff.
///
/// `weight01` is the rule topGO describes for its default algorithm: a term
/// is scored without the genes of descendants removed further down, then any
/// direct child at least as significant as the term gives up its genes too
/// and the term is rescored, until no remaining child is. (topGO 2.65's
/// `.sigRatio.01` counts every child as more significant, so topGO itself
/// removes every child's genes unconditionally.)
///
/// Removed genes leave the ranking for that term's test, as in topGO's elim,
/// rather than counting against it as non-members.
fn score_with_graph(
    scorer: &Scorer,
    hits: &[Option<Vec<u32>>],
    graph: &TermGraph,
    algorithm: Algorithm,
    options: &GseaOptions,
) -> Vec<Option<TermScore>> {
    let term_count = hits.len();
    let mut children = vec![Vec::<usize>::new(); term_count];
    for (child, parents) in graph.parents.iter().enumerate() {
        if hits[child].is_none() {
            continue;
        }
        for parent in parents {
            if hits[*parent].is_some() {
                children[*parent].push(child);
            }
        }
    }

    let mut order: Vec<usize> = (0..term_count)
        .filter(|term| hits[*term].is_some())
        .collect();
    order.sort_by_key(|term| std::cmp::Reverse(graph.ancestor_counts[*term]));

    let mut results: Vec<Option<GraphTerm>> = (0..term_count).map(|_| None).collect();
    for level in
        order.chunk_by(|left, right| graph.ancestor_counts[*left] == graph.ancestor_counts[*right])
    {
        let scored: Vec<(usize, GraphTerm)> = level
            .par_iter()
            .map(|term| {
                let term_hits = hits[*term].as_deref().expect("ordered terms are tested");
                let child_result =
                    |child: &usize| results[*child].as_ref().expect("children are scored first");
                let mut removed = Vec::new();
                for child in &children[*term] {
                    removed = union(&removed, &child_result(child).passed_up);
                }
                let result = match algorithm {
                    Algorithm::Elim => {
                        let score = scorer.score(*term, 0, term_hits, &removed);
                        let passed_up = if score.p_value <= options.elim_cutoff {
                            term_hits.to_vec()
                        } else {
                            removed
                        };
                        GraphTerm { score, passed_up }
                    }
                    Algorithm::Classic => unreachable!("classic terms skip the graph"),
                    Algorithm::Weight01 => {
                        let mut remaining = children[*term].clone();
                        let mut attempt = 0;
                        let mut score = scorer.score(*term, attempt, term_hits, &removed);
                        loop {
                            let (stronger, weaker): (Vec<usize>, Vec<usize>) =
                                remaining.into_iter().partition(|child| {
                                    child_result(child).score.at_least_as_significant(&score)
                                });
                            remaining = weaker;
                            if stronger.is_empty() {
                                break;
                            }
                            for child in &stronger {
                                let child_hits = hits[*child].as_deref().expect("child is tested");
                                removed = union(&removed, child_hits);
                            }
                            attempt += 1;
                            score = scorer.score(*term, attempt, term_hits, &removed);
                        }
                        GraphTerm {
                            score,
                            passed_up: removed,
                        }
                    }
                };
                (*term, result)
            })
            .collect();
        for (term, result) in scored {
            results[term] = Some(result);
        }
    }

    results
        .into_iter()
        .map(|result| result.map(|result| result.score))
        .collect()
}

fn union(left: &[u32], right: &[u32]) -> Vec<u32> {
    let mut merged = Vec::with_capacity(left.len() + right.len());
    let (mut i, mut j) = (0, 0);
    while i < left.len() && j < right.len() {
        match left[i].cmp(&right[j]) {
            std::cmp::Ordering::Less => {
                merged.push(left[i]);
                i += 1;
            }
            std::cmp::Ordering::Greater => {
                merged.push(right[j]);
                j += 1;
            }
            std::cmp::Ordering::Equal => {
                merged.push(left[i]);
                i += 1;
                j += 1;
            }
        }
    }
    merged.extend_from_slice(&left[i..]);
    merged.extend_from_slice(&right[j..]);
    merged
}

#[derive(Clone, Debug)]
struct TermScore {
    set_size: u32,
    es: f64,
    nes: f64,
    p_value: f64,
    /// The KS null the term was tested against, for the pooled FDR.
    null: NullRef,
    /// Ranks in the full ranking.
    leading_edge: Vec<u32>,
    peak_rank: u32,
}

#[derive(Clone, Debug)]
enum NullRef {
    /// Wilcoxon, or a term left with no members.
    None,
    /// The whole-ranking null for this set size.
    Shared(usize),
    /// A null drawn for this term's reduced ranking.
    Own(Arc<NullDistribution>),
}

impl TermScore {
    /// At least as significant as `other`: a smaller p-value, or the same
    /// p-value (as at the KS permutation floor) and at least as large |NES|.
    fn at_least_as_significant(&self, other: &Self) -> bool {
        self.p_value < other.p_value
            || (self.p_value == other.p_value && self.nes.abs() >= other.nes.abs())
    }
}

/// Per-ranking state shared by every term test.
struct Scorer<'a> {
    ranking: &'a Ranking,
    options: &'a GseaOptions,
    query_seed: u64,
    /// KS step weight at each rank.
    weights: Vec<f64>,
    /// Ascending-score midrank at each rank, over the whole ranking.
    midranks: Vec<f64>,
    /// Sum of `t^3 - t` over tie groups in the whole ranking.
    tie_total: f64,
    /// KS nulls over the whole ranking, by set size.
    nulls: BTreeMap<usize, NullDistribution>,
}

impl<'a> Scorer<'a> {
    fn new(
        query_index: usize,
        ranking: &'a Ranking,
        hits: &[Option<Vec<u32>>],
        options: &'a GseaOptions,
    ) -> Self {
        let query_seed = options.seed ^ (query_index as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15);
        let weights: Vec<f64> = ranking
            .scores
            .iter()
            .map(|score| step_weight(*score, options.weight))
            .collect();

        let n = ranking.len();
        let mut midranks = vec![0.0; n];
        let mut tie_total = 0.0;
        let mut start = 0;
        while start < n {
            let mut end = start + 1;
            while end < n && ranking.scores[end] == ranking.scores[start] {
                end += 1;
            }
            let ties = (end - start) as f64;
            // Ranks count up from the lowest score.
            let midrank = (n - end) as f64 + (ties + 1.0) / 2.0;
            midranks[start..end].fill(midrank);
            tie_total += ties * ties * ties - ties;
            start = end;
        }

        let nulls = if options.statistic == Statistic::Ks {
            let mut sizes: Vec<usize> = hits.iter().flatten().map(Vec::len).collect();
            sizes.sort_unstable();
            sizes.dedup();
            sizes
                .par_iter()
                .map(|k| {
                    let seed = query_seed ^ (*k as u64).wrapping_mul(0xD1B5_4A32_D192_ED03);
                    (
                        *k,
                        NullDistribution::sample(*k, &weights, options.permutations, seed),
                    )
                })
                .collect()
        } else {
            BTreeMap::new()
        };

        Self {
            ranking,
            options,
            query_seed,
            weights,
            midranks,
            tie_total,
            nulls,
        }
    }

    /// Scores one term with the ranks in `removed` (a sorted subset of `hits`)
    /// taken out of the ranking. `attempt` separates the random streams of
    /// repeated tests of one term.
    fn score(&self, target_index: usize, attempt: u32, hits: &[u32], removed: &[u32]) -> TermScore {
        let kept = difference(hits, removed);
        let n = self.ranking.len() - removed.len();
        let k = kept.len();
        if k == 0 || k >= n {
            return TermScore {
                set_size: k as u32,
                es: 0.0,
                nes: 0.0,
                p_value: 1.0,
                null: NullRef::None,
                leading_edge: Vec::new(),
                peak_rank: 0,
            };
        }
        // Ranks within the reduced ranking.
        let local = if removed.is_empty() {
            kept.clone()
        } else {
            reindex(&kept, removed)
        };

        let (es, nes, p_value, extreme, null) = match self.options.statistic {
            Statistic::Ks => {
                let (extreme, own) = if removed.is_empty() {
                    let extreme = running_sum_extreme(&local, Some(&self.weights), n);
                    (extreme, None)
                } else {
                    let weights = without(&self.weights, removed);
                    let extreme = running_sum_extreme(&local, Some(&weights), n);
                    let seed = self.query_seed
                        ^ splitmix((target_index as u64) << 16 | u64::from(attempt));
                    let null =
                        NullDistribution::sample(k, &weights, self.options.permutations, seed);
                    (extreme, Some(Arc::new(null)))
                };
                let (p_value, nes) = match own.as_deref() {
                    Some(null) => null.test(extreme.es),
                    None => self.nulls[&k].test(extreme.es),
                };
                let null = own.map_or(NullRef::Shared(k), NullRef::Own);
                (extreme.es, nes, p_value, extreme, null)
            }
            Statistic::Wilcoxon => {
                let (rank_sum, tie_total) = if removed.is_empty() {
                    (
                        kept.iter().map(|rank| self.midranks[*rank as usize]).sum(),
                        self.tie_total,
                    )
                } else {
                    self.reduced_rank_sum(&kept, removed)
                };
                let (es, z, p_value) = rank_sum_test(rank_sum, tie_total, k, n);
                let extreme = running_sum_extreme(&local, None, n);
                (es, z, p_value, extreme, NullRef::None)
            }
        };

        let (leading_edge, peak_rank) = if es >= 0.0 {
            (&kept[..=extreme.max_hit], kept[extreme.max_hit] + 1)
        } else {
            (&kept[extreme.min_hit..], kept[extreme.min_hit] + 1)
        };
        TermScore {
            set_size: k as u32,
            es,
            nes,
            p_value,
            null,
            leading_edge: leading_edge.to_vec(),
            peak_rank,
        }
    }

    /// Wilcoxon rank sum of `kept` and the tie total, over the ranking
    /// without `removed`.
    fn reduced_rank_sum(&self, kept: &[u32], removed: &[u32]) -> (f64, f64) {
        let scores = &self.ranking.scores;
        let n = scores.len();
        let reduced = n - removed.len();
        let (mut next_kept, mut next_removed) = (0, 0);
        let mut seen = 0;
        let mut rank_sum = 0.0;
        let mut tie_total = 0.0;
        let mut start = 0;
        while start < n {
            let mut end = start + 1;
            while end < n && scores[end] == scores[start] {
                end += 1;
            }
            let mut members = 0usize;
            while next_kept < kept.len() && (kept[next_kept] as usize) < end {
                members += 1;
                next_kept += 1;
            }
            let mut gone = 0usize;
            while next_removed < removed.len() && (removed[next_removed] as usize) < end {
                gone += 1;
                next_removed += 1;
            }
            let ties = end - start - gone;
            if ties > 0 {
                let ties_f = ties as f64;
                let midrank = (reduced - seen - ties) as f64 + (ties_f + 1.0) / 2.0;
                rank_sum += members as f64 * midrank;
                tie_total += ties_f * ties_f * ties_f - ties_f;
                seen += ties;
            }
            start = end;
        }
        (rank_sum, tie_total)
    }
}

/// Mann-Whitney U from the members' rank sum, as scipy's `mannwhitneyu`
/// computes it with `method="asymptotic"` and continuity correction. Returns
/// the rank-biserial correlation, the z-score, and the two-sided p-value.
fn rank_sum_test(rank_sum: f64, tie_total: f64, k: usize, n: usize) -> (f64, f64, f64) {
    let k_f = k as f64;
    let n_f = n as f64;
    let others = n_f - k_f;
    let u = rank_sum - k_f * (k_f + 1.0) / 2.0;
    let mean = k_f * others / 2.0;
    let variance = k_f * others / 12.0 * ((n_f + 1.0) - tie_total / (n_f * (n_f - 1.0)));
    let es = 2.0 * u / (k_f * others) - 1.0;
    if variance <= 0.0 {
        return (es, 0.0, 1.0);
    }
    let deviation = u - mean;
    let corrected = (deviation.abs() - 0.5).max(0.0).copysign(deviation);
    let z = corrected / variance.sqrt();
    let p_value = erfc(z.abs() / std::f64::consts::SQRT_2).min(1.0);
    (es, z, p_value)
}

/// Complementary error function, through the regularized upper incomplete
/// gamma function `Q(1/2, x^2)`; accurate to about 1e-15 relative across the
/// range where it does not underflow.
fn erfc(x: f64) -> f64 {
    if x < 0.0 {
        return 2.0 - erfc(-x);
    }
    const A: f64 = 0.5;
    // ln(Gamma(1/2)) = ln(sqrt(pi)).
    const LN_GAMMA_A: f64 = 0.572_364_942_924_700_1;
    let z = x * x;
    if z == 0.0 {
        return 1.0;
    }
    let prefactor = (-z + A * z.ln() - LN_GAMMA_A).exp();
    if z < A + 1.0 {
        // Series for the lower function P.
        let mut term = 1.0 / A;
        let mut sum = term;
        let mut ap = A;
        for _ in 0..500 {
            ap += 1.0;
            term *= z / ap;
            sum += term;
            if term.abs() < sum.abs() * f64::EPSILON {
                break;
            }
        }
        1.0 - sum * prefactor
    } else {
        // Lentz's continued fraction for the upper function Q.
        const TINY: f64 = 1e-300;
        let mut b = z + 1.0 - A;
        let mut c = 1.0 / TINY;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..500 {
            let i = f64::from(i);
            let an = -i * (i - A);
            b += 2.0;
            d = an * d + b;
            if d.abs() < TINY {
                d = TINY;
            }
            c = b + an / c;
            if c.abs() < TINY {
                c = TINY;
            }
            d = 1.0 / d;
            let delta = d * c;
            h *= delta;
            if (delta - 1.0).abs() < f64::EPSILON {
                break;
            }
        }
        prefactor * h
    }
}

/// Sorted `all` without the sorted subset `removed`.
fn difference(all: &[u32], removed: &[u32]) -> Vec<u32> {
    let mut next = 0;
    all.iter()
        .copied()
        .filter(|rank| {
            while next < removed.len() && removed[next] < *rank {
                next += 1;
            }
            !(next < removed.len() && removed[next] == *rank)
        })
        .collect()
}

/// Positions of the sorted `ranks` once the sorted, disjoint `removed` ranks
/// are taken out of the ranking.
fn reindex(ranks: &[u32], removed: &[u32]) -> Vec<u32> {
    let mut before = 0;
    ranks
        .iter()
        .map(|rank| {
            while before < removed.len() && removed[before] < *rank {
                before += 1;
            }
            rank - before as u32
        })
        .collect()
}

/// `values` without the entries at the sorted `removed` indices.
fn without(values: &[f64], removed: &[u32]) -> Vec<f64> {
    let mut next = 0;
    values
        .iter()
        .enumerate()
        .filter_map(|(index, value)| {
            if next < removed.len() && removed[next] as usize == index {
                next += 1;
                None
            } else {
                Some(*value)
            }
        })
        .collect()
}

fn splitmix(value: u64) -> u64 {
    let mut x = value.wrapping_add(0x9E37_79B9_7F4A_7C15);
    x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
    x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
    x ^ (x >> 31)
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
/// `None` weights every member equally.
fn running_sum_extreme(hits: &[u32], weights: Option<&[f64]>, n: usize) -> RunningSumExtreme {
    let k = hits.len();
    let hit_weight = |rank: u32| weights.map_or(1.0, |weights| weights[rank as usize]);
    let hit_total: f64 = hits.iter().map(|rank| hit_weight(*rank)).sum();
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
            hit_weight(*rank) / hit_total
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
#[derive(Debug)]
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
            let es = running_sum_extreme(&hits, Some(weights), n).es;
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

    /// Null NES magnitudes: each sign's ES over that sign's mean.
    fn normalized(&self) -> (Vec<f64>, Vec<f64>) {
        let positive_mean = mean(&self.positive);
        let negative_mean = mean(&self.negative).abs();
        (
            self.positive.iter().map(|es| es / positive_mean).collect(),
            self.negative
                .iter()
                .map(|es| es.abs() / negative_mean)
                .collect(),
        )
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

/// GSEA's FDR q-value (Subramanian et al. 2005): for an observed NES, the
/// share of the pooled null NES at least as extreme, over the share of
/// observed NES at least as extreme, within each sign. Pooling the normalized
/// nulls of every tested term resolves far below one term's permutation
/// floor. q-values are made monotone in |NES|, as BH is in p.
fn gsea_fdr(rows: &mut [GseaRow], nulls: &[NullRef], shared: &BTreeMap<usize, NullDistribution>) {
    let mut shared_uses = BTreeMap::<usize, f64>::new();
    for null in nulls {
        if let NullRef::Shared(k) = null {
            *shared_uses.entry(*k).or_default() += 1.0;
        }
    }
    let mut pooled: [Vec<(f64, f64)>; 2] = [Vec::new(), Vec::new()];
    let mut add = |null: &NullDistribution, uses: f64| {
        let (positive, negative) = null.normalized();
        pooled[0].extend(positive.into_iter().map(|value| (value, uses)));
        pooled[1].extend(negative.into_iter().map(|value| (value, uses)));
    };
    for (k, uses) in &shared_uses {
        add(&shared[k], *uses);
    }
    for null in nulls {
        if let NullRef::Own(null) = null {
            add(null, 1.0);
        }
    }

    for (sign, pool) in pooled.iter_mut().enumerate() {
        pool.sort_by(|left, right| left.0.total_cmp(&right.0));
        // Weight of null values at or above each position.
        let mut at_or_above = vec![0.0; pool.len() + 1];
        for index in (0..pool.len()).rev() {
            at_or_above[index] = at_or_above[index + 1] + pool[index].1;
        }
        let null_total = at_or_above[0];

        let mut members: Vec<usize> = (0..rows.len())
            .filter(|index| {
                let row = &rows[*index];
                !matches!(nulls[*index], NullRef::None)
                    && row.nes.is_finite()
                    && (row.es >= 0.0) == (sign == 0)
            })
            .collect();
        members.sort_by(|left, right| rows[*left].nes.abs().total_cmp(&rows[*right].nes.abs()));
        let observed_total = members.len() as f64;
        let mut first_tie = 0;
        let mut running: f64 = 1.0;
        for (position, index) in members.iter().enumerate() {
            let magnitude = rows[*index].nes.abs();
            if position > 0 && rows[members[position - 1]].nes.abs() < magnitude {
                first_tie = position;
            }
            let first = pool.partition_point(|(value, _)| *value < magnitude);
            let null_share = at_or_above[first] / null_total;
            let observed_share = (observed_total - first_tie as f64) / observed_total;
            let fdr = if null_total > 0.0 {
                null_share / observed_share
            } else {
                1.0
            };
            running = running.min(fdr.min(1.0));
            rows[*index].p_adjust = running;
        }
    }
    for (row, null) in rows.iter_mut().zip(nulls) {
        if matches!(null, NullRef::None) || !row.nes.is_finite() {
            row.p_adjust = 1.0;
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
        Algorithm, GseaOptions, RankedList, Ranking, Scorer, Statistic, benjamini_hochberg, erfc,
        gsea_one, running_sum_extreme,
    };
    use crate::index::{
        GeneUniverse, NamedSet, RawSet, TermGraph, build_flat_sets, build_ontology_sets,
        build_term_graph,
    };

    fn approx(left: f64, right: f64) -> bool {
        (left - right).abs() < 1e-12
    }

    fn close(left: f64, right: f64, relative: f64) -> bool {
        (left - right).abs() <= relative * right.abs()
    }

    #[test]
    fn running_sum_matches_hand_computation() {
        // Ten genes, scores 10..1, members at ranks 0, 2 and 7.
        let weights: Vec<f64> = (1..=10).rev().map(f64::from).collect();
        let hits = [0, 2, 7];
        let extreme = running_sum_extreme(&hits, Some(&weights), weights.len());
        // NR = 10 + 8 + 3; peak after rank 2: (10 + 8) / 21 - 1 / 7.
        assert!(approx(extreme.es, 18.0 / 21.0 - 1.0 / 7.0));
        assert_eq!(extreme.max_hit, 1);

        // Unweighted, members at the bottom: the trough is right before them.
        let extreme = running_sum_extreme(&[7, 8, 9], None, 10);
        assert!(approx(extreme.es, -1.0));
        assert_eq!(extreme.min_hit, 0);
    }

    #[test]
    fn erfc_matches_reference_values() {
        // scipy.special.erfc
        for (x, expected) in [
            (0.5, 0.479_500_122_186_953_5),
            (1.0, 0.157_299_207_050_285_16),
            (2.0, 0.004_677_734_981_047_266),
            (5.0, 1.537_459_794_428_034_7e-12),
            (10.0, 2.088_487_583_762_544_6e-45),
            (26.0, 5.663_192_408_856_145e-296),
        ] {
            assert!(close(erfc(x), expected, 1e-13), "erfc({x}) = {}", erfc(x));
        }
        assert!(approx(erfc(0.0), 1.0));
        assert!(close(erfc(-1.0), 2.0 - 0.157_299_207_050_285_16, 1e-15));
    }

    fn ranking_from(scores: &[f64]) -> (GeneUniverse, Ranking) {
        let genes: Vec<String> = (0..scores.len()).map(|i| format!("g{i:03}")).collect();
        let mut universe = GeneUniverse::new();
        universe.add_genes(&genes);
        let ranking = Ranking::new(&RankedList::new("r", genes, scores.to_vec()), &universe);
        (universe, ranking.unwrap())
    }

    #[test]
    fn wilcoxon_matches_scipy_with_ties_and_removal() {
        let scores = [
            5.0, 4.0, 4.0, 3.5, 3.0, 3.0, 3.0, 2.0, 1.5, 1.0, 0.5, 0.5, 0.0, -0.5, -1.0, -1.0,
            -2.0, -2.5, -3.0, -4.0,
        ];
        let (_, ranking) = ranking_from(&scores);
        let options = GseaOptions {
            statistic: Statistic::Wilcoxon,
            min_size: 1,
            ..GseaOptions::default()
        };
        let scorer = Scorer::new(0, &ranking, &[], &options);
        // scipy.stats.mannwhitneyu(method="asymptotic", use_continuity=True)
        let up = scorer.score(0, 0, &[0, 2, 4, 5, 9, 13], &[]);
        assert!(close(up.p_value, 0.057_168_546_991_619_5, 1e-12));
        assert!(approx(up.es, 0.559_523_809_523_809_5));
        let reduced = scorer.score(0, 0, &[0, 2, 4, 5, 9, 13], &[0, 2]);
        assert_eq!(reduced.set_size, 4);
        assert!(close(reduced.p_value, 0.311_525_077_078_600_33, 1e-12));
        assert!(approx(reduced.es, 0.357_142_857_142_857_2));
        let down = scorer.score(0, 0, &[14, 15, 17, 19], &[]);
        assert!(close(down.p_value, 0.012_051_584_623_132_423, 1e-12));
        assert!(approx(down.es, -0.843_75) && down.nes < 0.0);
    }

    fn fixture() -> (GeneUniverse, Ranking, Vec<NamedSet>) {
        let scores: Vec<f64> = (0..200).map(|i| 100.0 - i as f64).collect();
        let (universe, ranking) = ranking_from(&scores);
        let pick = |range: std::ops::Range<usize>, step: usize| -> Vec<String> {
            range.step_by(step).map(|i| format!("g{i:03}")).collect()
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
        for statistic in [Statistic::Ks, Statistic::Wilcoxon] {
            let options = GseaOptions {
                statistic,
                permutations: 2000,
                include_leading_edge: true,
                ..GseaOptions::default()
            };
            let rows = gsea_one(&ranking, &sets, None, &options).unwrap();
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
            assert_eq!(bottom.peak_rank, 161);

            assert!(row("spread").p_value > 0.1);
            assert_eq!(
                rows.last().map(|row| sets[row.target_index].id.as_str()),
                Some("spread")
            );
        }
    }

    #[test]
    fn is_reproducible_and_respects_size_limits() {
        let (_, ranking, sets) = fixture();
        let options = GseaOptions {
            permutations: 200,
            ..GseaOptions::default()
        };
        let first = gsea_one(&ranking, &sets, None, &options).unwrap();
        let second = gsea_one(&ranking, &sets, None, &options).unwrap();
        let summary =
            |rows: &[super::GseaRow]| rows.iter().map(|row| row.p_value).collect::<Vec<_>>();
        assert_eq!(summary(&first), summary(&second));

        let limited = GseaOptions {
            min_size: 21,
            ..options
        };
        assert!(
            gsea_one(&ranking, &sets, None, &limited)
                .unwrap()
                .is_empty()
        );
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
        let mut rows = gsea_one(&ranking, &sets, None, &GseaOptions::default()).unwrap();
        for (row, p_value) in rows.iter_mut().zip([0.01, 0.04, 0.03]) {
            row.p_value = p_value;
        }
        benjamini_hochberg(&mut rows);
        let adjusted: Vec<f64> = rows.iter().map(|row| row.p_adjust).collect();
        assert!(approx(adjusted[0], 0.03));
        assert!(approx(adjusted[1], 0.04));
        assert!(approx(adjusted[2], 0.04));
    }

    /// Two-term ontologies over a 200-gene ranking, scores falling with rank.
    fn ontology(
        annotations: &[(std::ops::Range<usize>, usize, &str)],
        closure: &[(&str, &str)],
    ) -> (Ranking, Vec<NamedSet>, TermGraph) {
        let scores: Vec<f64> = (0..200).map(|i| 100.0 - i as f64).collect();
        let (universe, ranking) = ranking_from(&scores);
        let annotations: Vec<(String, String)> = annotations
            .iter()
            .flat_map(|(range, step, term)| {
                range
                    .clone()
                    .step_by(*step)
                    .map(|i| (format!("g{i:03}"), (*term).to_owned()))
            })
            .collect();
        let closure: Vec<(String, String)> = closure
            .iter()
            .map(|(child, ancestor)| ((*child).to_owned(), (*ancestor).to_owned()))
            .collect();
        let sets = build_ontology_sets(&[], &closure, &annotations, &universe);
        let graph = build_term_graph(&sets, &closure);
        (ranking, sets, graph)
    }

    fn run(
        ranking: &Ranking,
        sets: &[NamedSet],
        graph: &TermGraph,
        statistic: Statistic,
        algorithm: Algorithm,
    ) -> Vec<super::GseaRow> {
        let options = GseaOptions {
            statistic,
            algorithm,
            min_size: 5,
            permutations: 2000,
            ..GseaOptions::default()
        };
        gsea_one(ranking, sets, Some(graph), &options).unwrap()
    }

    #[test]
    fn graph_algorithms_drop_a_parent_carried_by_its_child() {
        // The child holds 20 top genes; the parent adds 32 spread out.
        let (ranking, sets, graph) = ontology(
            &[(0..40, 2, "child"), (41..200, 5, "parent")],
            &[
                ("child", "child"),
                ("parent", "parent"),
                ("child", "parent"),
            ],
        );
        for statistic in [Statistic::Ks, Statistic::Wilcoxon] {
            let p_value = |algorithm, id: &str| {
                let rows = run(&ranking, &sets, &graph, statistic, algorithm);
                let row = rows.iter().find(|row| sets[row.target_index].id == id);
                let row = row.unwrap().clone();
                (row.p_value, row.eliminated, row.set_size)
            };
            assert!(p_value(Algorithm::Classic, "parent").0 < 0.01);
            for algorithm in [Algorithm::Elim, Algorithm::Weight01] {
                let (parent_p, eliminated, set_size) = p_value(algorithm, "parent");
                assert!(parent_p > 0.05, "{statistic:?} {algorithm:?} {parent_p}");
                assert_eq!((eliminated, set_size), (20, 32));
                assert!(p_value(algorithm, "child").0 < 0.01);
            }
        }
    }

    #[test]
    fn weight01_keeps_a_parent_stronger_than_its_children() {
        // Two children with interleaved, moderately high genes. Together they
        // are more significant than either, so weight01 keeps the parent while
        // elim, seeing significant children, empties it.
        let (ranking, sets, graph) = ontology(
            &[(0..120, 6, "left"), (3..123, 6, "right")],
            &[
                ("left", "left"),
                ("right", "right"),
                ("parent", "parent"),
                ("left", "parent"),
                ("right", "parent"),
            ],
        );
        let row = |algorithm, id: &str| {
            let rows = run(&ranking, &sets, &graph, Statistic::Wilcoxon, algorithm);
            rows.into_iter()
                .find(|row| sets[row.target_index].id == id)
                .unwrap()
        };
        let classic_parent = row(Algorithm::Classic, "parent");
        let left = row(Algorithm::Classic, "left");
        assert!(left.p_value < 0.01 && classic_parent.p_value < left.p_value);

        let weighted = row(Algorithm::Weight01, "parent");
        assert_eq!(weighted.eliminated, 0);
        assert_eq!(weighted.p_value, classic_parent.p_value);

        let elim = row(Algorithm::Elim, "parent");
        assert_eq!((elim.eliminated, elim.set_size, elim.p_value), (40, 0, 1.0));
    }
}
