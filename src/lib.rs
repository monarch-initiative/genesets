//! Fast, ontology-neutral gene set enrichment primitives.
//!
//! The crate is organized around a small set of data structures that are useful
//! both for ontology enrichment and for flat gene set libraries:
//!
//! - [`GeneUniverse`] assigns each gene a stable dense index.
//! - [`NamedSet`] stores a term, sample, or library entry as a dense bitset.
//! - [`enrichment`] scores query sets against target sets with Fisher exact
//!   enrichment.
//! - [`io`] loads simple TSV/CSV, GMT, and GMX input formats.
//!
//! The binary crate provides the supported CLI. The library surface is kept
//! small so benchmark and eval code can exercise the same core engine.

pub mod bitset;
pub mod enrichment;
pub mod fisher;
pub mod index;
pub mod io;
pub mod output;

pub use enrichment::{Correction, EnrichmentOptions, EnrichmentRow};
pub use index::{GeneUniverse, NamedSet, RawSet};
pub use io::SetFormat;
