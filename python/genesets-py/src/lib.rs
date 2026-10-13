//! Python bindings for genesets-rs.

use genesets_rs::{
    GeneUniverse, RawSet,
    index::{
        build_flat_sets, build_gene_list_bitset, build_ontology_sets, build_term_graph, union_sets,
    },
    mgsa::{BetaPrior, MgsaOptions, mgsa_one},
};
use pyo3::{exceptions::PyValueError, prelude::*, types::PyDict};

/// Runs MGSA for one study set and returns a dict with `terms` (highest
/// posterior first), `stats`, `study_size`, and `background_size`.
///
/// Results do not depend on the order of the inputs; `genes` lists are
/// sorted.
///
/// `annotations` are (gene, term) pairs that are already propagated to
/// ancestor terms. `closure` is optional (child, ancestor) pairs; it lets the
/// sampler propose moves between related terms and picks the most specific
/// member of a group of equivalent terms.
#[pyfunction]
#[pyo3(signature = (
    annotations,
    study,
    background = None,
    closure = None,
    *,
    samples_per_term = 100,
    burn_per_term = 10,
    term_prob = 0.5,
    term_count = 0.0,
    false_pos_prob = 0.5,
    false_pos_count = 0.0,
    false_neg_prob = 0.5,
    false_neg_count = 0.0,
    flip_rate = 1.0,
    step_rate = 1.0,
    jump_rate = 1.0,
    seed = 123_456_789,
))]
#[allow(clippy::too_many_arguments)]
fn mgsa<'py>(
    py: Python<'py>,
    annotations: Vec<(String, String)>,
    study: Vec<String>,
    background: Option<Vec<String>>,
    closure: Option<Vec<(String, String)>>,
    samples_per_term: u32,
    burn_per_term: u32,
    term_prob: f64,
    term_count: f64,
    false_pos_prob: f64,
    false_pos_count: f64,
    false_neg_prob: f64,
    false_neg_count: f64,
    flip_rate: f64,
    step_rate: f64,
    jump_rate: f64,
    seed: u64,
) -> PyResult<Bound<'py, PyDict>> {
    for prob in [term_prob, false_pos_prob, false_neg_prob] {
        if !(0.0..=1.0).contains(&prob) {
            return Err(PyValueError::new_err("prior modes must be between 0 and 1"));
        }
    }
    for count in [term_count, false_pos_count, false_neg_count] {
        if !count.is_finite() || count < 0.0 {
            return Err(PyValueError::new_err(
                "prior pseudocounts must be finite and non-negative",
            ));
        }
    }
    let rates = [flip_rate, step_rate, jump_rate];
    if rates.iter().any(|rate| *rate < 0.0) || rates.iter().all(|rate| *rate == 0.0) {
        return Err(PyValueError::new_err(
            "move rates must be non-negative and not all zero",
        ));
    }

    let options = MgsaOptions {
        term_prior: BetaPrior {
            mode: term_prob,
            count: term_count,
        },
        false_pos_prior: BetaPrior {
            mode: false_pos_prob,
            count: false_pos_count,
        },
        false_neg_prior: BetaPrior {
            mode: false_neg_prob,
            count: false_neg_count,
        },
        samples_per_term,
        burn_per_term,
        flip_rate,
        step_rate,
        jump_rate,
        seed,
        include_overlap_genes: true,
    };

    let (result, universe, targets) = py.allow_threads(|| {
        let mut universe = GeneUniverse::new();
        for (gene, _) in &annotations {
            universe.add_gene(gene);
        }
        universe.add_genes(&study);
        if let Some(background) = &background {
            universe.add_genes(background);
        }
        // Annotations arrive propagated, so no closure is applied here.
        let targets = build_ontology_sets(&[], &[], &annotations, &universe);
        let graph = closure
            .as_deref()
            .filter(|closure| !closure.is_empty())
            .map(|closure| build_term_graph(&targets, closure));
        let background_bits = background.as_deref().map_or_else(
            || union_sets(&targets, universe.len()),
            |genes| build_gene_list_bitset(genes, &universe),
        );
        let query = build_flat_sets(&[RawSet::new("study", None, study)], &universe).remove(0);
        let result = mgsa_one(&query, &targets, &background_bits, graph.as_ref(), &options);
        (result, universe, targets)
    });

    let terms = result
        .terms
        .iter()
        .map(|term| {
            let row = PyDict::new(py);
            row.set_item("id", &targets[term.target_index].id)?;
            row.set_item("posterior", term.posterior)?;
            row.set_item("overlap", term.overlap)?;
            row.set_item("term_size", term.target_size)?;
            row.set_item("p_value", term.p_value)?;
            row.set_item(
                "equivalent_ids",
                term.equivalent_indices
                    .iter()
                    .map(|index| targets[*index].id.as_str())
                    .collect::<Vec<_>>(),
            )?;
            let mut genes: Vec<&str> = term
                .overlap_indices
                .as_deref()
                .unwrap_or(&[])
                .iter()
                .map(|index| universe.gene_id(*index))
                .collect();
            genes.sort_unstable();
            row.set_item("genes", genes)?;
            Ok(row)
        })
        .collect::<PyResult<Vec<_>>>()?;

    let stats = &result.stats;
    let stats_dict = PyDict::new(py);
    stats_dict.set_item("modeled_terms", stats.modeled_terms)?;
    stats_dict.set_item("samples", stats.samples)?;
    stats_dict.set_item("burn", stats.burn)?;
    stats_dict.set_item(
        "acceptance_rate",
        if stats.samples > 0 {
            stats.accepted_moves as f64 / stats.samples as f64
        } else {
            0.0
        },
    )?;
    stats_dict.set_item("mean_active_terms", stats.mean_active_terms)?;
    stats_dict.set_item("term_rate", stats.term_rate)?;
    stats_dict.set_item("false_pos_rate", stats.false_pos_rate)?;
    stats_dict.set_item("false_neg_rate", stats.false_neg_rate)?;

    let out = PyDict::new(py);
    out.set_item("terms", terms)?;
    out.set_item("stats", stats_dict)?;
    out.set_item("study_size", result.query_size)?;
    out.set_item("background_size", result.background_size)?;
    Ok(out)
}

#[pymodule]
fn _genesets(module: &Bound<'_, PyModule>) -> PyResult<()> {
    module.add_function(wrap_pyfunction!(mgsa, module)?)?;
    module.add("__version__", env!("CARGO_PKG_VERSION"))?;
    Ok(())
}
