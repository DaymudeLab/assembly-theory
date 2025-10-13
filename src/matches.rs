//! Pairs of non-overlapping, isomorphic subgraphs in a molecule.

use std::collections::{BTreeMap, HashMap, HashSet};

use bit_set::BitSet;
use dashmap::DashMap;
use petgraph::graph::EdgeIndex;
use rayon::iter::{IntoParallelRefIterator, ParallelIterator};

use crate::{
    assembly::ParallelMode,
    canonize::{canonize, CanonizeMode, Labeling},
    enumerate::{enumerate_subgraphs, EnumerateMode},
    molecule::Molecule,
    state::State,
    utils::edge_neighbors,
};

struct DagNode {
    fragment: BitSet,
    canonical_id: usize,
    children: Vec<usize>,
}

/// Pairs of non-overlapping, isomorphic subgraphs in a molecule, sorted to
/// guarantee deterministic iteration.
pub struct Matches {
    dag: Vec<DagNode>,
    id_to_match: Vec<(usize, usize)>,
    match_to_id: HashMap<(usize, usize), usize>,
}

impl DagNode {
    pub fn new(fragment: BitSet, canonical_id: usize) -> Self {
        Self {
            fragment,
            canonical_id,
            children: Vec::new(),
        }
    }

    pub fn fragment(&self) -> &BitSet {
        &self.fragment
    }

    pub fn len(&self) -> usize {
        self.fragment.len()
    }
}

impl Matches {
    /// Generate [`Matches`] from the given molecule with the specified modes.
    pub fn new(
        mol: &Molecule,
        enumerate_mode: EnumerateMode,
        canonize_mode: CanonizeMode,
        parallel_mode: ParallelMode,
    ) -> Self {
        let num_edges = mol.graph().edge_count();

        let mut subgraph_ids: Vec<usize> = Vec::new();
        let mut constructed_frags: HashSet<BitSet> = HashSet::new();
        let mut dag_ids: HashMap<BitSet, usize> = HashMap::new();
        let mut dag: Vec<DagNode> = Vec::with_capacity(num_edges);
        let mut matches: Vec<(usize, usize)> = Vec::new();
        let mut next_canonical_id = 1;

        // Helper functions
        fn add_to_dag(
            dag: &mut Vec<DagNode>,
            fragment: &BitSet,
            parent_id: usize,
            canonical_id: usize,
        ) -> usize {
            // Create new dag node
            let new_dag_node = DagNode::new(fragment.clone(), canonical_id);
            let child_id = dag.len();

            // Add child to dag list
            dag.push(new_dag_node);

            // Add index to parent's child list
            let parent = &mut dag[parent_id];
            parent.children.push(child_id);

            child_id
        }

        fn create_buckets(
            mol: &Molecule,
            subgraph_ids: &Vec<usize>,
            dag: &Vec<DagNode>,
            constructed_frags: &mut HashSet<BitSet>,
            canonize_mode: CanonizeMode,
        ) -> BTreeMap<Labeling, Vec<(BitSet, usize)>> {
            let mut buckets: BTreeMap<Labeling, Vec<(BitSet, usize)>> = BTreeMap::new();
            let num_edges = mol.graph().edge_count();

            // Extend and bucket new subgraphs
            for parent_id in subgraph_ids.iter() {
                // Get fragment to be extended
                let fragment = &dag[*parent_id].fragment;

                // Construct edge neighborhood of fragment
                let mut neighborhood = BitSet::with_capacity(num_edges);
                for e in fragment {
                    neighborhood
                        .extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                }
                neighborhood.difference_with(fragment);

                // Extend fragment with edges from neighborhood
                for e in neighborhood.iter() {
                    let mut new_fragment = fragment.clone();
                    new_fragment.insert(e);

                    // Check if fragment has already been created
                    // If yes, continue
                    if constructed_frags.contains(&new_fragment) {
                        continue;
                    } else {
                        constructed_frags.insert(new_fragment.clone());
                    }

                    // Bucket subgraph
                    let canonical_id = canonize(mol, &new_fragment, canonize_mode);
                    buckets
                        .entry(canonical_id)
                        .and_modify(|bucket| bucket.push((new_fragment.clone(), *parent_id)))
                        .or_insert(vec![(new_fragment, *parent_id)]);
                }
            }

            buckets
        }

        // Generate subgraphs with one edge
        for i in 0..num_edges {
            let mut new_bitset = BitSet::with_capacity(num_edges);
            new_bitset.insert(i);

            let node = DagNode::new(new_bitset, 0);
            dag.push(node);
            subgraph_ids.push(i);
        }

        // Generate subgraphs with one more edge than the previous iteration
        while subgraph_ids.len() > 0 {
            let mut child_subgraph_ids: Vec<usize> = Vec::new();
            let buckets = create_buckets(
                mol,
                &subgraph_ids,
                &dag,
                &mut constructed_frags,
                canonize_mode,
            );

            // Go through buckets to create matches
            for isomorphic_frags in buckets.values() {
                // Variables to track if a match has been found
                let mut frag_has_match = BitSet::with_capacity(isomorphic_frags.len());

                // Loop through pairs of fragments
                for frag1_id in 0..isomorphic_frags.len() {
                    for frag2_id in frag1_id + 1..isomorphic_frags.len() {
                        let (frag1, parent1_id) = &isomorphic_frags[frag1_id];
                        let (frag2, parent2_id) = &isomorphic_frags[frag2_id];

                        if frag1.is_disjoint(&frag2) {
                            // If this is the first match for a fragment, create a dag node
                            if !frag_has_match.contains(frag1_id) {
                                let child_id =
                                    add_to_dag(&mut dag, &frag1, *parent1_id, next_canonical_id);
                                child_subgraph_ids.push(child_id);
                                dag_ids.insert(frag1.clone(), child_id);
                            }
                            if !frag_has_match.contains(frag2_id) {
                                let child_id =
                                    add_to_dag(&mut dag, &frag2, *parent2_id, next_canonical_id);
                                child_subgraph_ids.push(child_id);
                                dag_ids.insert(frag1.clone(), child_id);
                            }

                            // Get fragment ids and add to matches
                            let frag1_dag_id = dag_ids.get(&frag1).unwrap();
                            let frag2_dag_id = dag_ids.get(&frag2).unwrap();

                            if frag1_dag_id > frag2_dag_id {
                                matches.push((*frag1_dag_id, *frag2_dag_id));
                            } else {
                                matches.push((*frag2_dag_id, *frag1_dag_id));
                            }

                            // Mark that these fragments have matches
                            frag_has_match.insert(frag1_id);
                            frag_has_match.insert(frag2_id);
                        }
                    }
                }

                // If there was a match, increment canon_id counter
                if !frag_has_match.is_empty() {
                    next_canonical_id += 1;
                }
            }

            // Update to the new subgraphs
            subgraph_ids = child_subgraph_ids;
        }

        // Sort matches
        // Larger frag_ids get a smaller match_id
        // Thus fragments with many edges have small ids
        matches.sort();
        matches.reverse();

        // Give ids to matches
        let mut next_match_id = 0;
        let mut id_to_match = Vec::new();
        let mut match_to_id: HashMap<(usize, usize), usize> = HashMap::default();

        for m in matches {
            id_to_match.push(m);
            match_to_id.insert(m, next_match_id);
            next_match_id += 1;
        }

        Self {
            dag,
            id_to_match,
            match_to_id,
        }
    }

    /// Return the number of matches.
    pub fn len(&self) -> usize {
        self.id_to_match.len()
    }

    /// Return `true` if there are no matches.
    pub fn is_empty(&self) -> bool {
        self.id_to_match.is_empty()
    }

    /// Return all matches that are later than the last-removed match in the
    /// given assembly state.
    pub fn later_matches(&self, state: &State) -> &[(BitSet, BitSet)] {
        // TODO
    }
}
