//! Strucutral information on matches in the molecular graph - pairs of
//! non-overlapping, isomorphic subgraphs

use std::collections::{BTreeMap, HashMap, HashSet};

use bit_set::BitSet;
use petgraph::graph::EdgeIndex;

use crate::{
    canonize::{canonize, CanonizeMode, Labeling},
    molecule::Molecule,
    state::State,
    utils::edge_neighbors,
};

/// Holds data for nodes of the dag
struct DagNode {
    /// Connected subgraph of the molecule
    fragment: BitSet,
    /// IDs for isomorphism. Two DagNodes have the same canonical_id iff their
    /// fragments are ismorphic.
    canonical_id: usize,
    /// List of indices of this nodes children. If A is a child of B then
    /// A.fragment is B.fragment with an additional edge.
    children: Vec<usize>,
}

/// Structural information on the molecules matches - pairs of nonoverlapping
/// isomorphic subgraph.
pub struct Matches {
    /// Each node of the dag holds information on a duplicatable subgraph (i.e.
    /// subgraphs such that there exists a disjoint ismorphic copy of the
    /// subgraph in the molecule).
    dag: Vec<DagNode>,
    /// List of all possible matches stored as a pair of fragment IDs (the
    /// index of the fragment in the dag).
    id_to_match: Vec<(usize, usize)>,
    /// Given a match (as a pair of fragment IDs) returns its match ID, i.e.,
    /// its index in id_to_match
    match_to_id: HashMap<(usize, usize), usize>,
}

impl DagNode {
    /// Create new [`DagNode`] object.
    pub fn new(fragment: BitSet, canonical_id: usize) -> Self {
        Self {
            fragment,
            canonical_id,
            children: Vec::new(),
        }
    }
}

impl Matches {
    /// Generate [`Matches`] from the given molecule with the specified modes.
    pub fn new(mol: &Molecule, canonize_mode: CanonizeMode) -> Self {
        let num_edges = mol.graph().edge_count();

        let mut subgraph_ids: Vec<usize> = Vec::new();
        let mut constructed_frags: HashSet<BitSet> = HashSet::new();
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
            subgraph_ids: &[usize],
            dag: &[DagNode],
            constructed_frags: &mut HashSet<BitSet>,
            canonize_mode: CanonizeMode,
        ) -> BTreeMap<Labeling, Vec<(BitSet, usize)>> {
            let mut buckets: BTreeMap<Labeling, Vec<(BitSet, usize)>> = BTreeMap::new();
            let num_edges = mol.graph().edge_count();

            for parent_id in subgraph_ids.iter() {
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
        while !subgraph_ids.is_empty() {
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
                let mut index_to_id = HashMap::<usize, usize>::new();

                // Loop through pairs of fragments
                for frag1_index in 0..isomorphic_frags.len() {
                    for frag2_index in frag1_index + 1..isomorphic_frags.len() {
                        let (frag1, parent1_id) = &isomorphic_frags[frag1_index];
                        let (frag2, parent2_id) = &isomorphic_frags[frag2_index];

                        if frag1.is_disjoint(frag2) {
                            // If this is the first match for a fragment, create a dag node
                            if !frag_has_match.contains(frag1_index) {
                                let child_id =
                                    add_to_dag(&mut dag, frag1, *parent1_id, next_canonical_id);
                                child_subgraph_ids.push(child_id);
                                index_to_id.insert(frag1_index, child_id);
                            }
                            if !frag_has_match.contains(frag2_index) {
                                let child_id =
                                    add_to_dag(&mut dag, frag2, *parent2_id, next_canonical_id);
                                child_subgraph_ids.push(child_id);
                                index_to_id.insert(frag2_index, child_id);
                            }

                            // Get fragment ids and add to matches
                            let frag1_dag_id = index_to_id.get(&frag1_index).unwrap();
                            let frag2_dag_id = index_to_id.get(&frag2_index).unwrap();

                            if frag1_dag_id > frag2_dag_id {
                                matches.push((*frag1_dag_id, *frag2_dag_id));
                            } else {
                                matches.push((*frag2_dag_id, *frag1_dag_id));
                            }

                            // Mark that these fragments have matches
                            frag_has_match.insert(frag1_index);
                            frag_has_match.insert(frag2_index);
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
        let mut id_to_match = Vec::new();
        let mut match_to_id: HashMap<(usize, usize), usize> = HashMap::default();

        for (next_match_id, m) in matches.iter().enumerate() {
            id_to_match.push(*m);
            match_to_id.insert(*m, next_match_id);
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
    pub fn later_matches(&self, state: &State) -> Vec<usize> {
        let state_fragments = state.fragments();
        let last_removed = state.last_removed();

        let mut valid_matches: Vec<usize> = Vec::new();
        let mut subgraphs: Vec<(usize, usize)> = Vec::new();

        // Helper function
        fn create_buckets(
            subgraphs: &[(usize, usize)],
            fragments: &[BitSet],
            dag: &[DagNode],
        ) -> BTreeMap<usize, Vec<(usize, usize)>> {
            let mut buckets: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();

            // Extend each subgraph and bucket
            for (dag_id, state_id) in subgraphs {
                let state_fragment = &fragments[*state_id];
                let children_ids = &dag[*dag_id].children;

                for child_id in children_ids {
                    let child_node = &dag[*child_id];

                    // Check if this extension is valid
                    // i.e. the extended subgraph is contained in a fragment of state
                    let child_fragment = &child_node.fragment;
                    let valid = child_fragment.is_subset(state_fragment);

                    // Add child fragment to bucket
                    if valid {
                        let canonical_id = &child_node.canonical_id;
                        buckets
                            .entry(*canonical_id)
                            .and_modify(|bucket| bucket.push((*child_id, *state_id)))
                            .or_insert(vec![(*child_id, *state_id)]);
                    }
                }
            }

            buckets
        }

        // Create subgraphs of size 1
        for (state_id, fragment) in state_fragments.iter().enumerate() {
            for edge_idx in fragment.iter() {
                subgraphs.push((edge_idx, state_id));
            }
        }

        while !subgraphs.is_empty() {
            let mut new_subgraphs = Vec::new();
            let mut buckets = create_buckets(&subgraphs, state_fragments, &self.dag);

            for isomorphic_frags in buckets.values_mut() {
                // Variables to track if a match has been found
                let mut has_match = BitSet::with_capacity(isomorphic_frags.len());

                // Loop over pairs of fragments
                for i in 0..isomorphic_frags.len() {
                    for j in i + 1..isomorphic_frags.len() {
                        let mut frag1 = (i, isomorphic_frags[i]);
                        let mut frag2 = (j, isomorphic_frags[j]);

                        // Order fragments
                        if frag1.1 .0 < frag2.1 .0 {
                            std::mem::swap(&mut frag1, &mut frag2);
                        }

                        let frag1_bucket_id = frag1.0;
                        let frag1_dag_id = frag1.1 .0;
                        let frag1_state_id = frag1.1 .1;

                        let frag2_bucket_id = frag2.0;
                        let frag2_dag_id = frag2.1 .0;
                        let frag2_state_id = frag2.1 .1;

                        // Check for valid match
                        if let Some(match_id) = self.match_to_id.get(&(frag1_dag_id, frag2_dag_id))
                        {
                            if *match_id as isize > last_removed {
                                valid_matches.push(*match_id);

                                // If this is the first time seeing that this frag has a match, add it to the new_subgraphs list
                                if !has_match.contains(frag1_bucket_id) {
                                    new_subgraphs.push((frag1_dag_id, frag1_state_id));
                                    has_match.insert(frag1_bucket_id);
                                }
                                if !has_match.contains(frag2_bucket_id) {
                                    new_subgraphs.push((frag2_dag_id, frag2_state_id));
                                    has_match.insert(frag2_bucket_id);
                                }
                            }
                        }
                    }
                }
            }

            subgraphs = new_subgraphs;
        }

        valid_matches.sort();
        valid_matches
    }

    /// Takes a match ID and returns the two disjoint isomoprhic fragments
    /// composing this match.
    pub fn match_fragments(&self, match_id: usize) -> (&BitSet, &BitSet) {
        let (frag1_dag_id, frag2_dag_id) = self.id_to_match[match_id];
        let frag1 = &self.dag[frag1_dag_id].fragment;
        let frag2 = &self.dag[frag2_dag_id].fragment;

        (frag1, frag2)
    }
}
