//! Model-based gene set analysis (MGSA).
//!
//! Instead of testing each term independently, MGSA explains the whole query
//! set with a small number of "active" terms. A gene annotated to at least one
//! active term is expected to be in the query; other genes are not. Genes that
//! break this rule are false positives or false negatives. The output for each
//! term is its posterior probability of being active, so a parent term that is
//! only enriched because of a child term is explained away instead of being
//! reported alongside it.
//!
//! The model follows Bauer, Gagneur & Robinson (2010), "GOing Bayesian:
//! model-based gene set analysis of genome-scale data", in the collapsed form
//! used by Holmes & Mungall's wtfgenes (<https://github.com/evoldoers/wtfgenes>):
//!
//! - one shared Bernoulli rate for terms being active, one false-positive rate
//!   for genes outside every active term, and one false-negative rate for genes
//!   inside an active term;
//! - each rate has a Beta prior and is integrated out analytically, so the
//!   sampler only walks over term on/off states;
//! - the sampler is Metropolis-Hastings with wtfgenes' moves: `flip` toggles
//!   one term, `step` swaps an active term for an ontology neighbor (parent,
//!   child, or sibling), and `jump` swaps an active term for any other term.
//!
//! Like wtfgenes, only terms with at least one query gene are modeled, and
//! terms with identical gene sets in the background are merged into one
//! variable. The most specific term (the one with the most ancestors) speaks
//! for the group and the others are reported as equivalents.

use std::collections::HashMap;

use rayon::prelude::*;

use crate::{
    bitset::DenseBitSet,
    fisher::FisherCache,
    index::{NamedSet, TermGraph},
};

/// A Beta prior given as its mode and its number of pseudocounts, as in
/// wtfgenes. A count of zero gives the uniform Beta(1, 1).
#[derive(Copy, Clone, Debug, PartialEq)]
pub struct BetaPrior {
    pub mode: f64,
    pub count: f64,
}

impl BetaPrior {
    fn alpha_beta(self) -> (f64, f64) {
        (
            self.mode * self.count + 1.0,
            (1.0 - self.mode) * self.count + 1.0,
        )
    }
}

impl Default for BetaPrior {
    fn default() -> Self {
        Self {
            mode: 0.5,
            count: 0.0,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MgsaOptions {
    pub term_prior: BetaPrior,
    pub false_pos_prior: BetaPrior,
    pub false_neg_prior: BetaPrior,
    /// Recorded MCMC samples per modeled term.
    pub samples_per_term: u32,
    /// Burn-in samples per modeled term.
    pub burn_per_term: u32,
    pub flip_rate: f64,
    /// Ignored when no term graph is available.
    pub step_rate: f64,
    pub jump_rate: f64,
    pub seed: u64,
    pub include_overlap_genes: bool,
}

impl Default for MgsaOptions {
    fn default() -> Self {
        Self {
            term_prior: BetaPrior::default(),
            false_pos_prior: BetaPrior::default(),
            false_neg_prior: BetaPrior::default(),
            samples_per_term: 100,
            burn_per_term: 10,
            flip_rate: 1.0,
            step_rate: 1.0,
            jump_rate: 1.0,
            seed: 123_456_789,
            include_overlap_genes: false,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MgsaTerm {
    /// The target that represents this group of equivalent targets.
    pub target_index: usize,
    /// Other targets with exactly the same genes in the background.
    pub equivalent_indices: Vec<usize>,
    pub overlap: u32,
    pub target_size: u32,
    /// Posterior probability that the term is active.
    pub posterior: f64,
    /// Term-for-term hypergeometric right tail, for comparison.
    pub p_value: f64,
    pub overlap_indices: Option<Vec<usize>>,
}

#[derive(Clone, Debug, Default)]
pub struct MgsaStats {
    pub modeled_terms: usize,
    pub samples: u64,
    pub burn: u64,
    pub accepted_moves: u64,
    pub mean_active_terms: f64,
    /// Posterior means of the three integrated-out rates.
    pub term_rate: f64,
    pub false_pos_rate: f64,
    pub false_neg_rate: f64,
}

#[derive(Clone, Debug)]
pub struct MgsaResult {
    pub query_index: usize,
    pub query_size: u32,
    pub background_size: u32,
    /// Modeled terms, highest posterior first.
    pub terms: Vec<MgsaTerm>,
    pub stats: MgsaStats,
}

/// Runs MGSA independently for each query, in parallel.
pub fn mgsa_matrix(
    queries: &[NamedSet],
    targets: &[NamedSet],
    background: &DenseBitSet,
    graph: Option<&TermGraph>,
    options: &MgsaOptions,
) -> Vec<MgsaResult> {
    let fisher = FisherCache::new(background.count_ones() as usize);
    queries
        .par_iter()
        .enumerate()
        .map(|(query_index, query)| {
            run_query(
                query_index,
                query,
                targets,
                background,
                graph,
                options,
                &fisher,
            )
        })
        .collect()
}

pub fn mgsa_one(
    query: &NamedSet,
    targets: &[NamedSet],
    background: &DenseBitSet,
    graph: Option<&TermGraph>,
    options: &MgsaOptions,
) -> MgsaResult {
    let fisher = FisherCache::new(background.count_ones() as usize);
    run_query(0, query, targets, background, graph, options, &fisher)
}

fn run_query(
    query_index: usize,
    query: &NamedSet,
    targets: &[NamedSet],
    background: &DenseBitSet,
    graph: Option<&TermGraph>,
    options: &MgsaOptions,
    fisher: &FisherCache,
) -> MgsaResult {
    let query_bits = query.bits.and(background);
    let query_size = query_bits.count_ones();
    let background_size = background.count_ones();
    let groups = group_relevant_targets(&query_bits, targets, background, graph);
    let mut model = Model::new(&groups, &query_bits, background_size, graph, options);

    let modeled = groups.len() as u64;
    let samples = u64::from(options.samples_per_term) * modeled;
    let burn = u64::from(options.burn_per_term) * modeled;
    let step_rate = if graph.is_some() {
        options.step_rate
    } else {
        0.0
    };
    let rates = [options.flip_rate, step_rate, options.jump_rate];
    let mut rng = Rng::new(options.seed ^ (query_index as u64).wrapping_mul(0x9E37_79B9_7F4A_7C15));

    let mut occupancy = vec![0u64; groups.len()];
    let mut stats = MgsaStats {
        modeled_terms: groups.len(),
        samples,
        burn,
        ..MgsaStats::default()
    };
    if modeled > 0 && rates.iter().any(|rate| *rate > 0.0) {
        let mut log_likelihood = model.log_likelihood();
        for sample in 0..burn + samples {
            let accepted = model.sample_move(&rates, &mut rng, &mut log_likelihood);
            if sample < burn {
                continue;
            }
            stats.accepted_moves += u64::from(accepted);
            for term in &model.active_list {
                occupancy[*term as usize] += 1;
            }
            let (term_rate, false_pos_rate, false_neg_rate) = model.mean_rates();
            stats.mean_active_terms += model.active_list.len() as f64;
            stats.term_rate += term_rate;
            stats.false_pos_rate += false_pos_rate;
            stats.false_neg_rate += false_neg_rate;
        }
        if samples > 0 {
            let samples = samples as f64;
            stats.mean_active_terms /= samples;
            stats.term_rate /= samples;
            stats.false_pos_rate /= samples;
            stats.false_neg_rate /= samples;
        }
    }

    let mut terms: Vec<MgsaTerm> = groups
        .into_iter()
        .zip(occupancy)
        .map(|(group, occupied)| {
            let target = &targets[group.exemplar];
            let overlap = group.bits.and_count(&query_bits);
            let target_size = group.bits.count_ones();
            MgsaTerm {
                target_index: group.exemplar,
                equivalent_indices: group.equivalents,
                overlap,
                target_size,
                posterior: if samples > 0 {
                    occupied as f64 / samples as f64
                } else {
                    0.0
                },
                p_value: fisher.right_tail(background_size, target_size, query_size, overlap),
                overlap_indices: options
                    .include_overlap_genes
                    .then(|| query_bits.intersection_indices3(&target.bits, background)),
            }
        })
        .collect();
    terms.sort_by(|left, right| {
        right
            .posterior
            .total_cmp(&left.posterior)
            .then_with(|| left.p_value.total_cmp(&right.p_value))
            .then_with(|| {
                targets[left.target_index]
                    .id
                    .cmp(&targets[right.target_index].id)
            })
    });

    MgsaResult {
        query_index,
        query_size,
        background_size,
        terms,
        stats,
    }
}

/// Targets with identical genes in the background, modeled as one variable.
struct TermGroup {
    exemplar: usize,
    equivalents: Vec<usize>,
    members: Vec<usize>,
    bits: DenseBitSet,
}

fn group_relevant_targets(
    query_bits: &DenseBitSet,
    targets: &[NamedSet],
    background: &DenseBitSet,
    graph: Option<&TermGraph>,
) -> Vec<TermGroup> {
    let mut by_genes: HashMap<DenseBitSet, Vec<usize>> = HashMap::new();
    for (index, target) in targets.iter().enumerate() {
        if target.bits.and_count(query_bits) > 0 {
            by_genes
                .entry(target.bits.and(background))
                .or_default()
                .push(index);
        }
    }

    let specificity = |index: usize| graph.map_or(0, |graph| graph.ancestor_counts[index]);
    let mut groups: Vec<TermGroup> = by_genes
        .into_iter()
        .map(|(bits, mut members)| {
            members.sort_by(|left, right| {
                specificity(*right)
                    .cmp(&specificity(*left))
                    .then_with(|| targets[*left].id.cmp(&targets[*right].id))
            });
            let exemplar = members[0];
            let mut equivalents = members[1..].to_vec();
            equivalents.sort_by(|left, right| targets[*left].id.cmp(&targets[*right].id));
            TermGroup {
                exemplar,
                equivalents,
                members,
                bits,
            }
        })
        .collect();
    groups.sort_by_key(|group| group.exemplar);
    groups
}

/// Neighbors of each group for `step` moves: parents, children, and siblings
/// (other children of a parent). The relation is symmetric, which keeps the
/// Hastings ratio `|neighbors(from)| / |neighbors(to)|` correct.
fn group_neighbors(groups: &[TermGroup], graph: &TermGraph) -> Vec<Vec<u32>> {
    let mut group_of_target = HashMap::new();
    for (group_index, group) in groups.iter().enumerate() {
        for member in &group.members {
            group_of_target.insert(*member, group_index);
        }
    }

    let mut parents = vec![Vec::new(); groups.len()];
    let mut children = vec![Vec::new(); groups.len()];
    for (group_index, group) in groups.iter().enumerate() {
        for member in &group.members {
            for parent in &graph.parents[*member] {
                if let Some(&parent_group) = group_of_target.get(parent)
                    && parent_group != group_index
                {
                    parents[group_index].push(parent_group);
                    children[parent_group].push(group_index);
                }
            }
        }
    }
    for list in parents.iter_mut().chain(children.iter_mut()) {
        list.sort_unstable();
        list.dedup();
    }

    (0..groups.len())
        .map(|group_index| {
            let mut neighbors: Vec<u32> = parents[group_index]
                .iter()
                .chain(&children[group_index])
                .chain(
                    parents[group_index]
                        .iter()
                        .flat_map(|parent| &children[*parent]),
                )
                .filter(|other| **other != group_index)
                .map(|other| *other as u32)
                .collect();
            neighbors.sort_unstable();
            neighbors.dedup();
            neighbors
        })
        .collect()
}

/// Sufficient statistics of the collapsed likelihood.
#[derive(Copy, Clone, Debug, Default, PartialEq)]
struct Counts {
    terms_on: u64,
    terms_off: u64,
    /// Inactive genes in the query, and not in the query.
    false_pos: u64,
    true_neg: u64,
    /// Active genes not in the query, and in the query.
    false_neg: u64,
    true_pos: u64,
}

struct Model {
    term_genes: Vec<Vec<u32>>,
    neighbors: Vec<Vec<u32>>,
    in_query: Vec<bool>,
    active_terms_by_gene: Vec<u32>,
    active: Vec<bool>,
    active_list: Vec<u32>,
    active_position: Vec<u32>,
    counts: Counts,
    priors: [(f64, f64); 3],
}

impl Model {
    fn new(
        groups: &[TermGroup],
        query_bits: &DenseBitSet,
        background_size: u32,
        graph: Option<&TermGraph>,
        options: &MgsaOptions,
    ) -> Self {
        // Only genes in some modeled term can change state; the rest are
        // permanently inactive and enter the counts as a constant.
        let mut local_by_gene = HashMap::new();
        let mut in_query = Vec::new();
        let term_genes = groups
            .iter()
            .map(|group| {
                group
                    .bits
                    .ones()
                    .into_iter()
                    .map(|gene| {
                        *local_by_gene.entry(gene).or_insert_with(|| {
                            in_query.push(query_bits.contains(gene));
                            (in_query.len() - 1) as u32
                        })
                    })
                    .collect()
            })
            .collect();
        let neighbors = graph.map_or_else(
            || vec![Vec::new(); groups.len()],
            |graph| group_neighbors(groups, graph),
        );
        let query_size = u64::from(query_bits.count_ones());

        Self {
            term_genes,
            neighbors,
            active_terms_by_gene: vec![0; in_query.len()],
            in_query,
            active: vec![false; groups.len()],
            active_list: Vec::new(),
            active_position: vec![0; groups.len()],
            counts: Counts {
                terms_off: groups.len() as u64,
                false_pos: query_size,
                true_neg: u64::from(background_size) - query_size,
                ..Counts::default()
            },
            priors: [
                options.term_prior.alpha_beta(),
                options.false_pos_prior.alpha_beta(),
                options.false_neg_prior.alpha_beta(),
            ],
        }
    }

    /// Log marginal likelihood up to a constant, with the three rates
    /// integrated out against their Beta priors.
    fn log_likelihood(&self) -> f64 {
        let c = self.counts;
        let [term, false_pos, false_neg] = self.priors;
        ln_beta(term.0 + c.terms_on as f64, term.1 + c.terms_off as f64)
            + ln_beta(
                false_pos.0 + c.false_pos as f64,
                false_pos.1 + c.true_neg as f64,
            )
            + ln_beta(
                false_neg.0 + c.false_neg as f64,
                false_neg.1 + c.true_pos as f64,
            )
    }

    fn mean_rates(&self) -> (f64, f64, f64) {
        let c = self.counts;
        let [term, false_pos, false_neg] = self.priors;
        let mean = |(alpha, beta): (f64, f64), succ: u64, fail: u64| {
            (alpha + succ as f64) / (alpha + beta + (succ + fail) as f64)
        };
        (
            mean(term, c.terms_on, c.terms_off),
            mean(false_pos, c.false_pos, c.true_neg),
            mean(false_neg, c.false_neg, c.true_pos),
        )
    }

    fn toggle(&mut self, term: usize) {
        let on = !self.active[term];
        self.active[term] = on;
        if on {
            self.active_position[term] = self.active_list.len() as u32;
            self.active_list.push(term as u32);
            self.counts.terms_on += 1;
            self.counts.terms_off -= 1;
        } else {
            let position = self.active_position[term] as usize;
            self.active_list.swap_remove(position);
            if let Some(moved) = self.active_list.get(position) {
                self.active_position[*moved as usize] = position as u32;
            }
            self.counts.terms_on -= 1;
            self.counts.terms_off += 1;
        }

        let counts = &mut self.counts;
        for gene in &self.term_genes[term] {
            let gene = *gene as usize;
            let active_terms = &mut self.active_terms_by_gene[gene];
            let changed = if on {
                *active_terms += 1;
                *active_terms == 1
            } else {
                *active_terms -= 1;
                *active_terms == 0
            };
            if !changed {
                continue;
            }
            match (on, self.in_query[gene]) {
                (true, true) => {
                    counts.false_pos -= 1;
                    counts.true_pos += 1;
                }
                (true, false) => {
                    counts.true_neg -= 1;
                    counts.false_neg += 1;
                }
                (false, true) => {
                    counts.true_pos -= 1;
                    counts.false_pos += 1;
                }
                (false, false) => {
                    counts.false_neg -= 1;
                    counts.true_neg += 1;
                }
            }
        }
    }

    /// Proposes one move and accepts or rejects it. Returns whether the state
    /// changed.
    fn sample_move(&mut self, rates: &[f64; 3], rng: &mut Rng, log_likelihood: &mut f64) -> bool {
        let terms = self.active.len();
        let (toggles, log_hastings): (Vec<usize>, f64) = match rng.weighted_index(rates) {
            0 => (vec![rng.index(terms)], 0.0),
            kind => {
                if self.active_list.is_empty() {
                    return false;
                }
                let from = self.active_list[rng.index(self.active_list.len())] as usize;
                let to = if kind == 1 {
                    let neighbors = &self.neighbors[from];
                    if neighbors.is_empty() {
                        return false;
                    }
                    neighbors[rng.index(neighbors.len())] as usize
                } else {
                    rng.index(terms)
                };
                if self.active[to] {
                    return false;
                }
                let log_hastings = if kind == 1 {
                    (self.neighbors[from].len() as f64 / self.neighbors[to].len() as f64).ln()
                } else {
                    0.0
                };
                (vec![from, to], log_hastings)
            }
        };

        for term in &toggles {
            self.toggle(*term);
        }
        let proposed = self.log_likelihood();
        let log_ratio = proposed - *log_likelihood + log_hastings;
        if log_ratio >= 0.0 || rng.uniform() < log_ratio.exp() {
            *log_likelihood = proposed;
            true
        } else {
            for term in &toggles {
                self.toggle(*term);
            }
            false
        }
    }
}

fn ln_beta(alpha: f64, beta: f64) -> f64 {
    ln_gamma(alpha) + ln_gamma(beta) - ln_gamma(alpha + beta)
}

/// Lanczos approximation (g = 7, n = 9), accurate to about 1e-15 relative.
fn ln_gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if x < 0.5 {
        let pi = std::f64::consts::PI;
        return (pi / (pi * x).sin()).ln() - ln_gamma(1.0 - x);
    }
    let x = x - 1.0;
    let t = x + 7.5;
    let series = COEFFICIENTS[1..]
        .iter()
        .enumerate()
        .fold(COEFFICIENTS[0], |sum, (i, c)| {
            sum + c / (x + (i + 1) as f64)
        });
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + series.ln()
}

/// xoshiro256++ seeded through SplitMix64; small, fast, and reproducible.
struct Rng {
    state: [u64; 4],
}

impl Rng {
    fn new(seed: u64) -> Self {
        let mut z = seed;
        let mut next = || {
            z = z.wrapping_add(0x9E37_79B9_7F4A_7C15);
            let mut x = z;
            x = (x ^ (x >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
            x = (x ^ (x >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
            x ^ (x >> 31)
        };
        Self {
            state: [next(), next(), next(), next()],
        }
    }

    fn next_u64(&mut self) -> u64 {
        let s = &mut self.state;
        let result = s[0].wrapping_add(s[3]).rotate_left(23).wrapping_add(s[0]);
        let t = s[1] << 17;
        s[2] ^= s[0];
        s[3] ^= s[1];
        s[1] ^= s[2];
        s[0] ^= s[3];
        s[2] ^= t;
        s[3] = s[3].rotate_left(45);
        result
    }

    fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 * (1.0 / (1u64 << 53) as f64)
    }

    fn index(&mut self, len: usize) -> usize {
        ((u128::from(self.next_u64()) * len as u128) >> 64) as usize
    }

    fn weighted_index(&mut self, weights: &[f64]) -> usize {
        let total: f64 = weights.iter().sum();
        let mut remaining = self.uniform() * total;
        for (index, weight) in weights.iter().enumerate() {
            if *weight > 0.0 && remaining < *weight {
                return index;
            }
            remaining -= weight;
        }
        weights
            .iter()
            .rposition(|weight| *weight > 0.0)
            .unwrap_or(0)
    }
}

#[cfg(test)]
mod tests {
    use super::{BetaPrior, MgsaOptions, Model, group_relevant_targets, ln_gamma, mgsa_one};
    use crate::index::{
        GeneUniverse, build_flat_sets, build_ontology_sets, build_term_graph, union_sets,
    };

    fn pairs(items: &[(&str, &str)]) -> Vec<(String, String)> {
        items
            .iter()
            .map(|(left, right)| ((*left).to_owned(), (*right).to_owned()))
            .collect()
    }

    #[test]
    fn ln_gamma_matches_factorials() {
        assert!((ln_gamma(1.0)).abs() < 1e-12);
        assert!((ln_gamma(5.0) - 24f64.ln()).abs() < 1e-12);
        assert!((ln_gamma(101.0) - (1..=100).map(|n| (n as f64).ln()).sum::<f64>()).abs() < 1e-9);
        assert!((ln_gamma(0.5) - std::f64::consts::PI.sqrt().ln()).abs() < 1e-12);
    }

    /// The wtfgenes MCMC test fixture: a small ontology of animals queried
    /// with the three mutant genes.
    #[test]
    fn sampler_matches_exact_posterior() {
        let closure = pairs(&[
            ("animal", "animal"),
            ("mutant", "mutant"),
            ("arachnid", "arachnid"),
            ("arachnid", "animal"),
            ("mammal", "mammal"),
            ("mammal", "animal"),
            ("spider", "spider"),
            ("spider", "arachnid"),
            ("spider", "animal"),
            ("primate", "primate"),
            ("primate", "mammal"),
            ("primate", "animal"),
            ("human", "human"),
            ("human", "primate"),
            ("human", "mammal"),
            ("human", "animal"),
            ("gorilla", "gorilla"),
            ("gorilla", "primate"),
            ("gorilla", "mammal"),
            ("gorilla", "animal"),
            ("spiderhuman", "spiderhuman"),
            ("spiderhuman", "arachnid"),
            ("spiderhuman", "human"),
            ("spiderhuman", "mutant"),
            ("spiderhuman", "primate"),
            ("spiderhuman", "mammal"),
            ("spiderhuman", "animal"),
        ]);
        let annotations = pairs(&[
            ("peter-parker", "spiderhuman"),
            ("may-parker", "spiderhuman"),
            ("socrates", "human"),
            ("charlotte", "spider"),
            ("king-kong", "gorilla"),
            ("king-kong", "mutant"),
        ]);
        let mut universe = GeneUniverse::new();
        for (gene, _) in &annotations {
            universe.add_gene(gene);
        }
        let targets = build_ontology_sets(&[], &closure, &annotations, &universe);
        let graph = build_term_graph(&targets, &closure);
        let background = union_sets(&targets, universe.len());
        let query = build_flat_sets(
            &[crate::index::RawSet::new(
                "mutants",
                None,
                vec![
                    "peter-parker".into(),
                    "may-parker".into(),
                    "king-kong".into(),
                ],
            )],
            &universe,
        )
        .remove(0);

        // Same informative priors as the wtfgenes test.
        let options = MgsaOptions {
            term_prior: BetaPrior {
                mode: 0.1,
                count: 10.0,
            },
            false_pos_prior: BetaPrior {
                mode: 1.0 / 6.0,
                count: 6.0,
            },
            false_neg_prior: BetaPrior {
                mode: 1.0 / 6.0,
                count: 6.0,
            },
            samples_per_term: 40_000,
            burn_per_term: 100,
            ..MgsaOptions::default()
        };

        // Exact marginals by enumerating every term state.
        let query_bits = query.bits.and(&background);
        let groups = group_relevant_targets(&query_bits, &targets, &background, Some(&graph));
        let mut model = Model::new(
            &groups,
            &query_bits,
            background.count_ones(),
            Some(&graph),
            &options,
        );
        let n = groups.len();
        let mut log_likelihoods = Vec::with_capacity(1 << n);
        for state in 0u32..(1 << n) {
            for (term, active) in model.active.clone().into_iter().enumerate() {
                if active != (state & (1 << term) != 0) {
                    model.toggle(term);
                }
            }
            log_likelihoods.push(model.log_likelihood());
        }
        let max = log_likelihoods.iter().copied().fold(f64::MIN, f64::max);
        let weights: Vec<f64> = log_likelihoods.iter().map(|ll| (ll - max).exp()).collect();
        let norm: f64 = weights.iter().sum();
        let exact: Vec<f64> = (0..n)
            .map(|term| {
                weights
                    .iter()
                    .enumerate()
                    .filter(|(state, _)| state & (1 << term) != 0)
                    .map(|(_, weight)| weight / norm)
                    .sum()
            })
            .collect();

        let result = mgsa_one(&query, &targets, &background, Some(&graph), &options);
        assert_eq!(result.stats.modeled_terms, n);
        for (group, exact) in groups.iter().zip(&exact) {
            let term = result
                .terms
                .iter()
                .find(|term| term.target_index == group.exemplar)
                .unwrap();
            assert!(
                (term.posterior - exact).abs() < 0.02,
                "{}: sampled {} exact {}",
                targets[group.exemplar].id,
                term.posterior,
                exact
            );
        }

        // "mutant" alone covers every query gene, so its ancestors and
        // descendants are explained away.
        let top = &targets[result.terms[0].target_index].id;
        assert_eq!(top, "mutant");
        let animal = result
            .terms
            .iter()
            .find(|term| targets[term.target_index].id == "animal")
            .unwrap();
        assert!(animal.posterior < 0.1);
    }

    #[test]
    fn merges_equivalent_terms_into_most_specific() {
        let closure = pairs(&[
            ("child", "child"),
            ("child", "parent"),
            ("parent", "parent"),
        ]);
        let annotations = pairs(&[("g1", "child"), ("g2", "child")]);
        let mut universe = GeneUniverse::new();
        for gene in ["g1", "g2", "g3"] {
            universe.add_gene(gene);
        }
        let targets = build_ontology_sets(&[], &closure, &annotations, &universe);
        let graph = build_term_graph(&targets, &closure);
        let background = union_sets(&targets, universe.len());
        let query = build_flat_sets(
            &[crate::index::RawSet::new("q", None, vec!["g1".into()])],
            &universe,
        )
        .remove(0);

        let result = mgsa_one(
            &query,
            &targets,
            &background,
            Some(&graph),
            &MgsaOptions::default(),
        );
        assert_eq!(result.terms.len(), 1);
        assert_eq!(targets[result.terms[0].target_index].id, "child");
        assert_eq!(
            result.terms[0]
                .equivalent_indices
                .iter()
                .map(|index| targets[*index].id.as_str())
                .collect::<Vec<_>>(),
            vec!["parent"]
        );
    }
}
