//! Strucutral information on matches in a molecular graph.

use std::collections::{BTreeMap, HashMap, HashSet};

use bit_set::BitSet;
use petgraph::graph::EdgeIndex;

use crate::{
    bounds::{match_bounds, Bound},
    canonize::{canonize, CanonizeMode, Labeling},
    molecule::Molecule,
    state::State,
    utils::{connected_components_under_edges, edge_neighbors},
};

/// A node in the DAG storing fragment information; see [`Matches`].
struct DagNode {
    /// The fragment (i.e., connected molecular subgraph) this node represents.
    fragment: BitSet,
    /// The canonical ID of this node's fragment. Two [`DagNode`]s have the
    /// same canonical ID iff their fragments are isomorphic.
    canonical_id: usize,
    /// Indices of this node's children/out-neighbors in the DAG. If u is a
    /// child of v then u.fragment is v.fragment with an additional edge.
    children: Vec<usize>,
}

/// Structural information on "matches" in a molecular graph, i.e., pairs of
/// edge-disjoint, isomorphic subgraphs.
pub struct Matches {
    /// Seet et al. (2024) perform match enumeration by constructing a directed
    /// acyclic graph (DAG). Each node in this DAG is a fragment (i.e., a
    /// connected molecular subgraph) that is duplicatable (i.e., there exists
    /// some other edge-disjoint fragment it is isomorphic to). If there is an
    /// edge from u to v, then fragment v is fragment u with one added edge
    /// (note: this new edge may be between two existing nodes in u or "extend"
    /// from one existing node in u to a new node).
    dag: Vec<DagNode>,
    /// All possible matches (i.e., pairs of edge-disjoint, isomorphic
    /// fragments) stored as pairs of fragment (i.e., DAG node) indices.
    matches: Vec<(usize, usize)>,
    /// Maps pairs of matched fragment indices to their indices in `matches`.
    match_to_ix: HashMap<(usize, usize), usize>,
}

impl DagNode {
    /// Create a new [`DagNode`].
    pub fn new(fragment: BitSet, canonical_id: usize) -> Self {
        Self {
            fragment,
            canonical_id,
            children: Vec::new(),
        }
    }
}

impl Matches {
    /// Generate [`Matches`] from the given molecule and canonization mode.
    pub fn new(mol: &Molecule, canonize_mode: CanonizeMode) -> Self {
        let num_edges = mol.graph().edge_count();
        let mut dag: Vec<DagNode> = Vec::with_capacity(num_edges);
        let mut matches: Vec<(usize, usize)> = Vec::new();

        // Generate all singleton edge fragments and add them to the DAG.
        let mut parent_frag_ixs: Vec<usize> = Vec::new();
        for i in 0..num_edges {
            // Create a fragment for the i-th edge in the molecule.
            let mut frag = BitSet::with_capacity(num_edges);
            frag.insert(i);

            // Add the fragment to the DAG. Since all singleton edge fragments
            // are trivially isomorphic, give them the same canonical ID.
            dag.push(DagNode::new(frag, 0));
            parent_frag_ixs.push(i);
        }

        // Iteratively extend all parent fragments by a single edge to form
        // child fragments. Add any child fragment that has an edge-disjoint
        // isomorphic match to the DAG.
        let mut explored_frags: HashSet<BitSet> = HashSet::new();
        let mut next_canonical_id = 1;
        while !parent_frag_ixs.is_empty() {
            // Extend all parent fragments by one edge and bin them into
            // isomorphism classes using canonization.
            let mut isomorphism_classes: BTreeMap<Labeling, Vec<(BitSet, usize)>> = BTreeMap::new();
            for parent_frag_ix in parent_frag_ixs {
                // Get this parent fragment's indcident edges.
                let parent_frag = &dag[parent_frag_ix].fragment;
                let mut incident_edges = BitSet::with_capacity(num_edges);
                for e in parent_frag {
                    incident_edges
                        .extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                }
                incident_edges.difference_with(parent_frag);

                // Form child fragments by extending this parent fragment by
                // each of its incident edges one by one, avoiding redundancy.
                for e in &incident_edges {
                    let mut child_fragment = parent_frag.clone();
                    child_fragment.insert(e);

                    if explored_frags.insert(child_fragment.clone()) {
                        isomorphism_classes
                            .entry(canonize(mol, &child_fragment, canonize_mode))
                            .and_modify(|c| c.push((child_fragment.clone(), parent_frag_ix)))
                            .or_insert(vec![(child_fragment, parent_frag_ix)]);
                    }
                }
            }

            // Iterate through each isomorphism class, adding fragments forming
            // a match (i.e., an edge-disjoint isomorphic pair) to the DAG.
            let mut child_frag_ixs: Vec<usize> = Vec::new();
            for isomorphic_frags in isomorphism_classes.values() {
                // Track which fragments in this isomorphism class are matched
                // and what their fragment (DAG node) indices are.
                let mut frag_has_match = BitSet::with_capacity(isomorphic_frags.len());
                let mut iso_to_frag_ix = HashMap::<usize, usize>::new();

                // Check if pairs of isomorphic fragments form a match.
                for iso_ix1 in 0..isomorphic_frags.len() {
                    for iso_ix2 in (iso_ix1 + 1)..isomorphic_frags.len() {
                        let (frag1, frag1_parent_ix) = &isomorphic_frags[iso_ix1];
                        let (frag2, frag2_parent_ix) = &isomorphic_frags[iso_ix2];

                        // If fragments are edge-disjoint, they are a match.
                        if frag1.is_disjoint(frag2) {
                            // Add matched fragments to the DAG on first match.
                            if frag_has_match.insert(iso_ix1) {
                                let frag1_ix = dag.len();
                                dag.push(DagNode::new(frag1.clone(), next_canonical_id));
                                dag[*frag1_parent_ix].children.push(frag1_ix);
                                iso_to_frag_ix.insert(iso_ix1, frag1_ix);
                                child_frag_ixs.push(frag1_ix);
                            }
                            if frag_has_match.insert(iso_ix2) {
                                let frag2_ix = dag.len();
                                dag.push(DagNode::new(frag2.clone(), next_canonical_id));
                                dag[*frag2_parent_ix].children.push(frag2_ix);
                                iso_to_frag_ix.insert(iso_ix2, frag2_ix);
                                child_frag_ixs.push(frag2_ix);
                            }

                            // Store the match using the fragment indices.
                            let frag1_ix = iso_to_frag_ix.get(&iso_ix1).unwrap();
                            let frag2_ix = iso_to_frag_ix.get(&iso_ix2).unwrap();
                            if frag1_ix > frag2_ix {
                                matches.push((*frag1_ix, *frag2_ix));
                            } else {
                                matches.push((*frag2_ix, *frag1_ix));
                            }
                        }
                    }
                }

                // If there was a match, increment to the next canonical ID.
                if !frag_has_match.is_empty() {
                    next_canonical_id += 1;
                }
            }

            // Use the child fragments as the parents in the next iteration.
            parent_frag_ixs = child_frag_ixs;
        }

        // Sort matches in descending order of their fragment index pairs. This
        // puts matches containing fragments with the most edges first.
        matches.sort_by_key(|m| std::cmp::Reverse(*m));

        // Store matches and their indices for easy lookup later.
        let mut match_to_ix: HashMap<(usize, usize), usize> = HashMap::new();
        match_to_ix.extend(matches.iter().enumerate().map(|(ix, m)| (*m, ix)));

        Self {
            dag,
            matches,
            match_to_ix,
        }
    }

    /// Return the number of matches.
    pub fn len(&self) -> usize {
        self.matches.len()
    }

    /// Return `true` if there are no matches.
    pub fn is_empty(&self) -> bool {
        self.matches.is_empty()
    }

    /// Return all matches whose removal from the given assembly state may
    /// result in a better assembly index according to the given match bounds.
    ///
    /// A match (i.e., two edge-disjoint isomorphic fragments) is removable
    /// from an assembly state if (1) each match fragment is a subgraph of some
    /// assembly state fragment, and (2) the match's index is strictly greater
    /// than that of the last match removed from this assembly state.
    pub fn matches_to_remove(
        &self,
        mol: &Molecule,
        state: &State,
        best: usize,
        bounds: &[Bound],
    ) -> (Vec<BitSet>, Vec<usize>) {
        // The search for removable matches uses the DAG of duplicatable
        // fragments, starting with singleton edge fragments (the DAG's source
        // nodes). Because removable fragments are subgraphs of assembly state
        // fragments, we track both a fragment's index (in the DAG) and the
        // indices of assembly state fragments that contain it.
        let mut frag_state_ixs: Vec<(usize, usize)> = Vec::new();
        for (state_ix, frag) in state.fragments().iter().enumerate() {
            frag_state_ixs.extend(frag.iter().map(|e| (e, state_ix)));
        }

        // For later use when applying match bounds, create a container for
        // "matchable edge masks" indexed by match size (i.e., the number of
        // edges in each of a match's fragments). Edge masks for a given match
        // size indicate, for each fragment in this assembly state, which edges
        // are included in some match of the given size.
        let mut matchable_edge_masks: Vec<Vec<BitSet>> = Vec::new();
        let mut iso_classes_by_len: Vec<BTreeMap<usize, Vec<(usize, usize)>>> = Vec::new();

        // Collect duplicatable fragments by isomorphism and find which edges
        // are used in a removable match. Use BFS over the DAG, only extending
        // search paths that include fragments contained in removable matches.
        // (If some duplicatable fragment is not in a removable match, then no
        // child fragment that extends it with additional edges can be either.)
        while !frag_state_ixs.is_empty() {
            // Extend the search by one additional level in the DAG, binning
            // relevant child fragments into isomorphism classes.
            let mut isomorphism_classes: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
            for (frag_ix, state_ix) in frag_state_ixs {
                // Get this fragment's children in the DAG and the assembly
                // state's fragment that contains it.
                let child_frag_ixs = &self.dag[frag_ix].children;
                let state_frag = &state.fragments()[state_ix];

                // Collect all child fragments that are subgraphs of the state
                // fragment and bin them into isomorphism classes.
                for child_frag_ix in child_frag_ixs {
                    let child_dag_node = &self.dag[*child_frag_ix];
                    if child_dag_node.fragment.is_subset(state_frag) {
                        isomorphism_classes
                            .entry(child_dag_node.canonical_id)
                            .and_modify(|c| c.push((*child_frag_ix, state_ix)))
                            .or_insert(vec![(*child_frag_ix, state_ix)]);
                    }
                }
            }

            // Iterate through each isomorphism class (containing isomorphic
            // removable fragments) to identify removable matches.
            let mut next_frag_state_ixs: Vec<(usize, usize)> = Vec::new();
            let mut matchable_edges =
                vec![BitSet::with_capacity(mol.graph().edge_count()); state.fragments().len()];
            for isomorphic_frags in isomorphism_classes.values_mut() {
                // Track which fragments in this isomorphism class are matched.
                let mut frag_has_match = BitSet::with_capacity(isomorphic_frags.len());

                // Check if pairs of isomorphic fragments form a match.
                for iso_ix1 in 0..isomorphic_frags.len() {
                    for iso_ix2 in (iso_ix1 + 1)..isomorphic_frags.len() {
                        // If both fragments have already been found to be
                        // in some match(es), don't check this again.
                        if frag_has_match.contains(iso_ix1) && frag_has_match.contains(iso_ix2) {
                            continue;
                        }

                        // Order the fragments by descending fragment index.
                        let mut frag1_ixs = (isomorphic_frags[iso_ix1], iso_ix1);
                        let mut frag2_ixs = (isomorphic_frags[iso_ix2], iso_ix2);
                        if isomorphic_frags[iso_ix1].0 < isomorphic_frags[iso_ix2].0 {
                            std::mem::swap(&mut frag1_ixs, &mut frag2_ixs);
                        }

                        // Unpack fragment indices.
                        let ((frag1_ix, frag1_state_ix), frag1_iso_ix) = frag1_ixs;
                        let ((frag2_ix, frag2_state_ix), frag2_iso_ix) = frag2_ixs;

                        // If these fragments match and occur strictly later
                        // in the match order than the last match removed by
                        // this assembly state, they are a removable match.
                        if let Some(match_ix) = self.match_to_ix.get(&(frag1_ix, frag2_ix)) {
                            if *match_ix as isize > state.last_removed() {
                                // Extend the search for removable matches with
                                // these fragments on first match.
                                if frag_has_match.insert(frag1_iso_ix) {
                                    next_frag_state_ixs.push((frag1_ix, frag1_state_ix));
                                    matchable_edges[frag1_state_ix]
                                        .union_with(&self.dag[frag1_ix].fragment);
                                }
                                if frag_has_match.insert(frag2_iso_ix) {
                                    next_frag_state_ixs.push((frag2_ix, frag2_state_ix));
                                    matchable_edges[frag2_state_ix]
                                        .union_with(&self.dag[frag2_ix].fragment);
                                }
                            }
                        }
                    }
                }

                // Remove any fragments that were not part of a match.
                for ix in (0..isomorphic_frags.len())
                    .filter(|x| !frag_has_match.contains(*x))
                    .rev()
                {
                    isomorphic_frags.remove(ix);
                }
            }

            iso_classes_by_len.push(isomorphism_classes);
            matchable_edge_masks.push(matchable_edges);

            // Use the updated fragment/state indices in the next iteration.
            frag_state_ixs = next_frag_state_ixs;
        }

        // Breaking out of the loop implies that there were no matches found in
        // the last iteration. Thus the final entry of matchable_edge_masks
        // will be empty and can be discarded.
        matchable_edge_masks.pop();

        // Create new fragments for the current assembly state by removing
        // non-matched edges. Such edges will remain non-matchable lower in the
        // search tree and thus can be discarded.
        let intermediate_frags = {
            if !matchable_edge_masks.is_empty() {
                matchable_edge_masks[0]
                    .iter()
                    .flat_map(|frag| connected_components_under_edges(mol.graph(), frag))
                    .filter(|frag| frag.len() >= 2)
                    .collect::<Vec<BitSet>>()
            } else {
                vec![]
            }
        };

        // Get a lower bound on the size of the largest match whose removal may
        // yield a better assembly index.
        let match_bound = match_bounds(state.index(), best, &matchable_edge_masks, bounds);

        // Use the stored isomorphism classes to generate removable matches.
        let mut removable_matches: Vec<usize> = Vec::new();
        for (bucket_ix, isomorphism_classes) in iso_classes_by_len.iter().enumerate().rev() {
            // Stop if matches of this size are bounded by the match bounds.
            let match_len = bucket_ix + 2;
            if match_len < match_bound {
                break;
            }

            // Check if pairs of isomorphic fragments form a match.
            for isomorphic_frags in isomorphism_classes.values() {
                for i in 0..isomorphic_frags.len() {
                    for j in i + 1..isomorphic_frags.len() {
                        let mut frag1_id = isomorphic_frags[i].0;
                        let mut frag2_id = isomorphic_frags[j].0;

                        if frag1_id < frag2_id {
                            std::mem::swap(&mut frag1_id, &mut frag2_id);
                        }

                        if let Some(match_ix) = self.match_to_ix.get(&(frag1_id, frag2_id)) {
                            if *match_ix as isize > state.last_removed() {
                                removable_matches.push(*match_ix);
                            }
                        }
                    }
                }
            }
        }

        // Sort removable matches in ascending order of match index (i.e.,
        // those with larger fragments first).
        removable_matches.sort();
        (intermediate_frags, removable_matches)
    }

    /// Return the two edge-disjoint isomorphic fragments composing this match.
    pub fn match_fragments(&self, match_ix: usize) -> (&BitSet, &BitSet) {
        let (frag1_ix, frag2_ix) = self.matches[match_ix];

        (&self.dag[frag1_ix].fragment, &self.dag[frag2_ix].fragment)
    }
}
