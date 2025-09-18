use std::collections::{BTreeMap, HashMap, HashSet};

use bit_set::BitSet;
use fxhash::FxHashMap;
use itertools::Itertools;
use petgraph::graph::EdgeIndex;

use crate::{bounds::Bound, canonize::{canonize, CanonizeMode, Labeling}, molecule::{Bond, Element, Molecule}, reductions::CompatGraph, utils::{connected_components_under_edges, edge_neighbors}};

pub struct DagNode {
    fragment: BitSet,
    canon_id: usize,
    children: Vec<usize>,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Hash)]
struct EdgeType {
    bond: Bond,
    ends: (Element, Element),
}

pub struct Matches {
    id_to_match: Vec<(usize, usize)>,
    match_to_id: FxHashMap<(usize, usize), usize>,
    dag: Vec<DagNode>,
    // clique: CompatGraph,
    // clique_offset: usize,
    edge_types: Vec<usize>,
}


impl DagNode {
    pub fn new(fragment: BitSet, canon_id: usize) -> Self {
        Self {
            fragment,
            canon_id,
            children: Vec::new(),
        }
    }
}

impl Matches {
    // Generate DAG, clique, and other info from the base molecule.
    pub fn new(mol: &Molecule, canonize_mode: CanonizeMode) -> Self {
        let num_edges = mol.graph().edge_count();

        let mut subgraphs: Vec<usize> = Vec::new();
        let mut constructed_frags: HashSet<BitSet> = HashSet::new();
        let mut frag_to_id: HashMap<BitSet, usize> = HashMap::new();

        let mut dag: Vec<DagNode> = Vec::with_capacity(num_edges);
        let mut matches: Vec<(usize, usize)> = Vec::new();

        let mut next_canon_id = 0;

        // Generate subgraphs with one edge
        for i in 0..num_edges {
            let mut new_bitset = BitSet::with_capacity(num_edges);
            new_bitset.insert(i);
            
            let node = DagNode::new(new_bitset, 0);

            dag.push(node);
            subgraphs.push(i);
        }

        next_canon_id += 1;

        while subgraphs.len() > 0 {
            // Initialize new buckets and subgraph list
            let mut buckets: BTreeMap<Labeling, Vec<(BitSet, usize)>> = BTreeMap::new();
            let mut child_subgraphs: Vec<usize> = Vec::new();

            // Extend and bucket new subgraphs
            for parent_idx in subgraphs.iter() {
                // Get fragment to be extended
                let fragment = &dag[*parent_idx].fragment;

                // Construct edge neighborhood of fragment
                let mut neighborhood = BitSet::with_capacity(num_edges);
                for e in fragment {
                    neighborhood.extend(edge_neighbors(mol.graph(), EdgeIndex::new(e)).map(|x| x.index()));
                }
                neighborhood.difference_with(fragment);

                // Extend fragment with edges from neighborhood
                for e in neighborhood.iter() {
                    // Create extended fragment
                    let mut new_fragment = fragment.clone();
                    new_fragment.insert(e);

                    // Check if fragment has already been created
                    // If yes, continue
                    if constructed_frags.contains(&new_fragment) {
                        continue;
                    }

                    // Bucket subgraph
                    // Gives ownership of new_fragment to buckets
                    let canon = canonize(mol, &new_fragment, canonize_mode);
                    buckets.entry(canon)
                        .and_modify(|bucket| bucket.push((new_fragment.clone(), *parent_idx)))
                        .or_insert(vec![(new_fragment.clone(), *parent_idx)]);
                
                    // Add to constructed_frags
                    constructed_frags.insert(new_fragment);
                }
            }

            // Go through buckets to create matches
            for fragments in buckets.values() {
                // Variables to track if a match has been found
                let mut frag_has_match = BitSet::with_capacity(fragments.len());

                // Loop through pairs of bitsets to find matches
                for pair in fragments.iter().enumerate().combinations(2) {
                    // Extract pair of fragments
                    let frag1 = &pair[0].1.0;
                    let frag2 = &pair[1].1.0;

                    // Check for disjointness
                    // If the fragments are disjoint, they form a match
                    if frag1.is_disjoint(&frag2) {
                        // Extract indices of fragments
                        let idx1 = pair[0].0;
                        let idx2 = pair[1].0;

                        // If this is the first match for a fragment, create a dag node
                        if !frag_has_match.contains(idx1) {
                            // Create new dag node
                            let new_dag_node = DagNode::new(frag1.clone(), next_canon_id);
                            let child_idx = dag.len();

                            // Add child to dag list
                            dag.push(new_dag_node);

                            // Add index to parent's child list
                            let parent_idx = pair[0].1.1;
                            let parent = &mut dag[parent_idx];
                            parent.children.push(child_idx);

                            // Add dag node to next subgraphs list
                            child_subgraphs.push(child_idx);

                            // Add to frag_to_id map
                            frag_to_id.insert(frag1.clone(), child_idx);
                        }
                        if !frag_has_match.contains(idx2) {
                            // Create new dag node
                            let new_dag_node = DagNode::new(frag2.clone(), next_canon_id);
                            let child_idx = dag.len();

                            // Add child to dag list
                            dag.push(new_dag_node);

                            // Add index to parent's child list
                            let parent_idx = pair[1].1.1;
                            let parent = &mut dag[parent_idx];
                            parent.children.push(child_idx);

                            // Add dag node to next subgraphs list
                            child_subgraphs.push(child_idx);

                            // Add to frag_to_id map
                            frag_to_id.insert(frag2.clone(), child_idx);
                        }

                        // Get fragment ids and add to matches
                        let frag1_id = frag_to_id.get(&frag1).unwrap();
                        let frag2_id = frag_to_id.get(&frag2).unwrap();

                        if frag1_id > frag2_id {
                            matches.push((*frag1_id, *frag2_id));
                        }
                        else {
                            matches.push((*frag2_id, *frag1_id));
                        }

                        // Mark that these fragments have matches
                        frag_has_match.insert(idx1);
                        frag_has_match.insert(idx2);
                    }
                }

                // If there was a match, increment canon_id counter
                if !frag_has_match.is_empty() {
                    next_canon_id += 1;
                }
            }

            // Update to the new subgraphs
            subgraphs = child_subgraphs;
        }

        // Sort matches
        // Larger frag_ids get a smaller match_id
        // Thus fragments with many edges have small ids
        matches.sort();
        matches.reverse();

        // Give ids to matches
        let mut next_match_id = 0;
        let mut id_to_match = Vec::new();
        let mut match_to_id: FxHashMap<(usize, usize), usize> = FxHashMap::default();

        for m in matches {
            id_to_match.push(m);
            match_to_id.insert(m, next_match_id);
            next_match_id += 1;
        }

        // Create edge type ids
        let graph = mol.graph();
        let mut edgetype_to_id: HashMap<EdgeType, usize> = HashMap::new();
        let mut edge_types: Vec<usize> = vec![0; num_edges];
        let mut next_id = 0;

        // Create a list of the elements of the nodes
        // The list can be indexed by NodeIndex-es
        let mut node_elements: Vec<Element> = Vec::new();
        for v in graph.node_weights() {
            node_elements.push(v.element());
        }

        // Extract bond types
        let weights: Vec<Bond> = graph.edge_weights().copied().collect();

        // types will hold an element for every unique edge type in fragment
        for e in graph.edge_indices() {
            // Extract bond type and endpoint atom types
            let bond = weights[e.index()];
            let (end1, end2) = graph.edge_endpoints(e).expect("Edge Endpoint Error");

            let atom1 = node_elements[end1.index()];
            let atom2 = node_elements[end2.index()];

            // Create edgetype
            let ends = if atom1 < atom2 {(atom1, atom2)} else {(atom2, atom1)};
            let edgetype = EdgeType{bond, ends};

            // Find or create id for this edgetype
            // and store the edgetype id for this edge
            if let Some(id) = edgetype_to_id.get(&edgetype) {
                edge_types[e.index()] = *id;
            }
            else {
                edge_types[e.index()] = next_id;
                edgetype_to_id.insert(edgetype, next_id);
                next_id += 1;
            }
        }

        Self {
            id_to_match,
            match_to_id,
            dag,
            edge_types,
        }
    }

    pub fn generate_matches(&self, mol: &Molecule, state: &[BitSet], state_index: usize, last_removed: usize, best: usize, bounds: &[Bound]) -> Vec<(usize, usize)> {
        let num_edges = mol.graph().edge_count();

        let mut valid_matches: Vec<(usize, usize)> = Vec::new();
        let mut subgraphs: Vec<(usize, usize)> = Vec::new();
        let mut masks: Vec<Vec<BitSet>> = Vec::new();
        let mut buckets_by_len: Vec<HashMap<usize, Vec<(usize, usize)>>> = Vec::new();

        // Create subgraphs of size 1
        for (state_id, fragment) in state.iter().enumerate() {
            for edge_idx in fragment.iter() {
                subgraphs.push((edge_idx, state_id));
            }
        }

        while subgraphs.len() > 0 {
            // Subgraphs to extend in next loop
            let mut buckets: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
            let mut new_subgraphs = Vec::new();
            let mut used_edge_mask = vec![BitSet::with_capacity(num_edges); state.len()];

            // Extend each subgraph and bucket
            for (frag_id, state_id) in subgraphs {
                // Get state frament that this subgraph is contained in
                let state_frag = &state[state_id];

                // Get ids of extensions of this subgraph
                let children_ids = &self.dag[frag_id].children;

                for child_id in children_ids {
                    // Check if this extension is valid
                    // i.e. the extended subgraph is contained in a fragment of state
                    let child_frag = &self.dag[*child_id].fragment;
                    let valid = child_frag.is_subset(state_frag);

                    if valid {
                        let canon_id = &self.dag[*child_id].canon_id;
                        
                        // Add fragment to bucket
                        buckets.entry(*canon_id)
                            .and_modify(|bucket| bucket.push((*child_id, state_id)))
                            .or_insert(vec![(*child_id, state_id)]);
                    }
                }
            }

            // Search through buckets and create matches
            for fragments in buckets.values_mut() {
                let mut has_match = BitSet::with_capacity(fragments.len());
                // Loop over pairs of fragments
                for i in 0..fragments.len() {
                    for j in i+1..fragments.len() {
                        if has_match.contains(i) && has_match.contains(j) {
                            continue;
                        }

                        let mut frag1 = (i, fragments[i]);
                        let mut frag2 = (j, fragments[j]);

                        // Order fragments
                        if frag1.1.0 < frag2.1.0 {
                            let temp = frag1;
                            frag1 = frag2;
                            frag2 = temp;
                        }

                        let frag1_bucket_idx = frag1.0;
                        let frag1_id = frag1.1.0;
                        let frag1_state_id = frag1.1.1;

                        let frag2_bucket_idx = frag2.0;
                        let frag2_id = frag2.1.0;
                        let frag2_state_id = frag2.1.1;

                        // Check for valid match
                        if let Some(x) = self.match_to_id.get(&(frag1_id, frag2_id)) {
                            // Only add match if it occurs later in the ordering than
                            // the last removed match
                            if *x >= last_removed {
                                // valid_matches.push((frag1_id, frag2_id));

                                // If this is the first time seeing that this frag has a match,
                                // add it to the new_subgraphs list and modify mask
                                if !has_match.contains(frag1_bucket_idx) {
                                    new_subgraphs.push((frag1_id, frag1_state_id));
                                    used_edge_mask[frag1_state_id].union_with(&self.dag[frag1_id].fragment);
                                    has_match.insert(frag1_bucket_idx);
                                }
                                if !has_match.contains(frag2_bucket_idx) {
                                    new_subgraphs.push((frag2_id, frag2_state_id));
                                    used_edge_mask[frag2_state_id].union_with(&self.dag[frag2_id].fragment);
                                    has_match.insert(frag2_bucket_idx);
                                }
                            }
                        }
                    }
                }

                // Remove any matchless fragments from bucket
                let len = fragments.len();
                for i in (0..len).rev() {
                    if !has_match.contains(i) {
                        fragments.remove(i);
                    }
                }     
            }

            // Add buckets to list
            buckets_by_len.push(buckets);

            // Add mask to list
            masks.push(used_edge_mask);

            // Update to new subgraphs
            subgraphs = new_subgraphs;
        }

        // Let new state be only the edges which can be matched
        let mut state = Vec::new();

        for frag in masks[0].iter() {
            // Add connected components to new state
            state.extend(connected_components_under_edges(mol.graph(), &frag));
        }

        // Remove frags of size 1 or less
        state.retain(|frag| frag.len() >= 2);


        // Use bounding strategy to find the largest length match
        // to try removing from this state.
        let mut largest_length = 2;

        for (i, list) in masks.iter().enumerate() {
            let mut stop = false;
            for bound_type in bounds {
                if stop {
                    break;
                }

                match bound_type {
                    Bound::Int => {
                        let mut bound = 0;
                        let i = i + 2;
                        for (j, frag) in list.iter().enumerate() {
                            let mf = masks[0][j].len();
                            let mi = frag.len();
                            let x = mf + (mi % i) - mi;
                            bound += mf - (mi / i) - (x+i-2) / (i-1);
                        }
                        bound -= (i as f32).log2().ceil() as usize;

                        if state_index - bound >= best {
                            stop = true;
                        }
                    }
                    Bound::VecSimple => {
                        let size_two_list = &masks[0];
                        let mut total_set = HashSet::new();
                        let mut small_set = HashSet::new();
                        
                        for frag in size_two_list.iter() {
                            for edge in frag {
                                total_set.insert(self.edge_types[edge]);
                            }
                        }
                        for frag in list {
                            for edge in frag {
                                small_set.insert(self.edge_types[edge]);
                            }
                        }

                        let z = total_set.len();
                        let zi = small_set.len();
                        let i = i + 2;
                        let s: usize = size_two_list.iter().map(|frag| frag.len()).sum();
                        let mi: usize = list.iter().map(|frag| frag.len()).sum();
                        let mf = s - mi;

                        let bound = s - {
                            if z < i {
                                z-1 + (i as f32 / z as f32).log2().ceil() as usize + (mi/i) + ((mi%i)+mf - (z - zi) + i-2)/(i-1)
                            }
                            else {
                                z + (mi-zi)/i + ((mi-zi)%i + mf - (z - zi) + i-2)/(i-1)
                            }
                        };

                        if mi == 0 {
                            break;
                        }

                        if state_index - bound >= best {
                            stop = true;
                        }                    
                    }
                    _ => ()
                }
            }

            if stop {
                largest_length += 1;
            }
            else {
                break;
            }
        }

        // Create matches
        for bucket in buckets_by_len[largest_length-2..].iter().rev() {
            for fragments in bucket.values() {
                for i in 0..fragments.len() {
                    for j in i+1..fragments.len() {
                        let mut frag1_id = fragments[i].0;
                        let mut frag2_id = fragments[j].0;

                        if frag1_id < frag2_id {
                            let temp = frag1_id;
                            frag1_id = frag2_id;
                            frag2_id = temp;
                        }

                        if let Some(x) = self.match_to_id.get(&(frag1_id, frag2_id)){
                            if *x >= last_removed {
                                valid_matches.push((frag1_id, frag2_id));
                            }
                        }
                    }
                }
            }
        }

        // Sort matches
        valid_matches.sort();
        valid_matches.reverse();

        valid_matches
    }

    pub fn len(&self) -> usize {
        self.id_to_match.len()
    }

    pub fn get_frag(&self, id: usize) -> &BitSet {
        &self.dag[id].fragment
    }

    pub fn get_match_id(&self, mat: &(usize, usize)) -> Option<usize> {
        self.match_to_id.get(mat).copied()
    }
}

