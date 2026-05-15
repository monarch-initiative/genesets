use std::collections::HashMap;

use crate::bitset::DenseBitSet;

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct RawSet {
    pub id: String,
    pub name: Option<String>,
    pub genes: Vec<String>,
}

impl RawSet {
    pub fn new(id: impl Into<String>, name: Option<String>, genes: Vec<String>) -> Self {
        Self {
            id: id.into(),
            name,
            genes,
        }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NamedSet {
    pub id: String,
    pub name: Option<String>,
    pub bits: DenseBitSet,
}

#[derive(Clone, Debug, Default)]
pub struct GeneUniverse {
    ids: Vec<String>,
    names: Vec<Option<String>>,
    map: HashMap<String, usize>,
}

impl GeneUniverse {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_gene(&mut self, gene: &str) -> usize {
        if let Some(index) = self.map.get(gene) {
            return *index;
        }
        let index = self.ids.len();
        self.ids.push(gene.to_owned());
        self.names.push(None);
        self.map.insert(gene.to_owned(), index);
        index
    }

    pub fn add_genes<'a>(&mut self, genes: impl IntoIterator<Item = &'a String>) {
        for gene in genes {
            self.add_gene(gene);
        }
    }

    pub fn get(&self, gene: &str) -> Option<usize> {
        self.map.get(gene).copied()
    }

    pub fn len(&self) -> usize {
        self.ids.len()
    }

    pub fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    pub fn gene_id(&self, index: usize) -> &str {
        &self.ids[index]
    }

    pub fn gene_name(&self, index: usize) -> Option<&str> {
        self.names[index].as_deref()
    }

    pub fn apply_gene_names(&mut self, names: &HashMap<String, String>) {
        for (index, gene_id) in self.ids.iter().enumerate() {
            self.names[index] = names.get(gene_id).cloned();
        }
    }
}

#[derive(Default)]
struct TermIndex {
    ids: Vec<String>,
    names: Vec<Option<String>>,
    map: HashMap<String, usize>,
}

impl TermIndex {
    fn add_term(&mut self, term: &str) -> usize {
        if let Some(index) = self.map.get(term) {
            return *index;
        }
        let index = self.ids.len();
        self.ids.push(term.to_owned());
        self.names.push(None);
        self.map.insert(term.to_owned(), index);
        index
    }

    fn set_name(&mut self, term: &str, name: &str) {
        let index = self.add_term(term);
        self.names[index] = Some(name.to_owned());
    }

    fn get(&self, term: &str) -> Option<usize> {
        self.map.get(term).copied()
    }

    fn len(&self) -> usize {
        self.ids.len()
    }
}

pub fn build_ontology_sets(
    term_names: &[(String, String)],
    closure: &[(String, String)],
    annotations: &[(String, String)],
    universe: &GeneUniverse,
) -> Vec<NamedSet> {
    let mut terms = TermIndex::default();

    for (term, name) in term_names {
        terms.set_name(term, name);
    }
    for (_, term) in annotations {
        terms.add_term(term);
    }
    for (child, ancestor) in closure {
        terms.add_term(child);
        terms.add_term(ancestor);
    }

    let mut ancestors_by_child = vec![Vec::<usize>::new(); terms.len()];
    for (child, ancestor) in closure {
        let child_index = terms.get(child).expect("closure child term was indexed");
        let ancestor_index = terms
            .get(ancestor)
            .expect("closure ancestor term was indexed");
        ancestors_by_child[child_index].push(ancestor_index);
    }
    for (term_index, ancestors) in ancestors_by_child.iter_mut().enumerate() {
        ancestors.push(term_index);
        ancestors.sort_unstable();
        ancestors.dedup();
    }

    let mut sets: Vec<NamedSet> = (0..terms.len())
        .map(|index| NamedSet {
            id: terms.ids[index].clone(),
            name: terms.names[index].clone(),
            bits: DenseBitSet::new(universe.len()),
        })
        .collect();

    for (gene, term) in annotations {
        let Some(gene_index) = universe.get(gene) else {
            continue;
        };
        let Some(term_index) = terms.get(term) else {
            continue;
        };
        for ancestor_index in &ancestors_by_child[term_index] {
            sets[*ancestor_index].bits.set(gene_index);
        }
    }

    sets
}

pub fn build_flat_sets(raw_sets: &[RawSet], universe: &GeneUniverse) -> Vec<NamedSet> {
    raw_sets
        .iter()
        .map(|raw| {
            let mut bits = DenseBitSet::new(universe.len());
            for gene in &raw.genes {
                if let Some(index) = universe.get(gene) {
                    bits.set(index);
                }
            }
            NamedSet {
                id: raw.id.clone(),
                name: raw.name.clone(),
                bits,
            }
        })
        .collect()
}

pub fn build_gene_list_bitset(genes: &[String], universe: &GeneUniverse) -> DenseBitSet {
    let mut bits = DenseBitSet::new(universe.len());
    for gene in genes {
        if let Some(index) = universe.get(gene) {
            bits.set(index);
        }
    }
    bits
}

pub fn union_sets(sets: &[NamedSet], gene_count: usize) -> DenseBitSet {
    let mut union = DenseBitSet::new(gene_count);
    for set in sets {
        union.or_with(&set.bits);
    }
    union
}

#[cfg(test)]
mod tests {
    use super::{GeneUniverse, build_ontology_sets};

    #[test]
    fn propagates_annotations_to_ancestors() {
        let mut universe = GeneUniverse::new();
        for gene in ["g1", "g2", "g3"] {
            universe.add_gene(gene);
        }

        let annotations = vec![
            ("g1".to_owned(), "child".to_owned()),
            ("g2".to_owned(), "parent".to_owned()),
        ];
        let closure = vec![
            ("child".to_owned(), "child".to_owned()),
            ("child".to_owned(), "parent".to_owned()),
            ("parent".to_owned(), "parent".to_owned()),
        ];
        let sets = build_ontology_sets(&[], &closure, &annotations, &universe);
        let child = sets.iter().find(|set| set.id == "child").unwrap();
        let parent = sets.iter().find(|set| set.id == "parent").unwrap();

        assert_eq!(child.bits.count_ones(), 1);
        assert_eq!(parent.bits.count_ones(), 2);
    }
}
