//! Fast, ontology-neutral gene set enrichment primitives.
//!
//! The crate is organized around a small set of data structures that are useful
//! both for ontology enrichment and for flat gene set libraries:
//!
//! - [`GeneUniverse`] assigns each gene a stable dense index.
//! - [`NamedSet`] stores a term, sample, or library entry as a dense bitset.
//! - [`enrichment`] scores query sets against target sets with Fisher exact
//!   enrichment.
//! - [`mgsa`] scores terms jointly with Bayesian model-based gene set
//!   analysis, so redundant ancestors of an enriched term are explained away.
//! - [`gsea`] scores target sets against a ranked, scored gene list with
//!   preranked gene set enrichment analysis.
//! - [`io`] loads simple TSV/CSV, GMT, and GMX input formats, and ranked lists.
//!
//! The binary crate provides the supported CLI. The library surface is kept
//! small so benchmark and eval code can exercise the same core engine.

pub mod bitset;
pub mod compare;
pub mod enrichment;
pub mod fisher;
pub mod gsea;
pub mod index;
pub mod io;
pub mod mgsa;
pub mod output;
mod rng;

pub use enrichment::{Correction, EnrichmentOptions, EnrichmentRow};
pub use index::{GeneUniverse, NamedSet, RawSet};
pub use io::SetFormat;
