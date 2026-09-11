# This script depends on having already run ring-data.csv
# It will inspect the output of that script and generate a control dataset
# This dataset will have the same number of molecules, but all molecules will be free
# from rings entirely. Each molecule in this dataset will have a number of 
# bonds and double bonds corresponding to the number of bonds/double bonds in the ring data. 

import random
import networkx as nx
import pandas as pd
from rdkit import Chem
from rdkit.Chem import RWMol
from pathlib import Path

def random_tree_prufer(n, max_degree=4):
    """Generate a random tree with n nodes using Prüfer sequence. Reject if any node exceeds max_degree."""
    while True:
        seq = [random.randint(0, n - 1) for _ in range(n - 2)]
        deg = [1] * n
        for v in seq:
            deg[v] += 1
        if any(d > max_degree for d in deg):
            # print("rejecting sequence")
            continue
        
        G = nx.from_prufer_sequence(seq)
        
        return G


def assign_double_bonds(G, D, max_valence=4):
    """Randomly assign D double bonds to edges of G while obeying valence constraint."""
    edges = list(G.edges())
    if D > len(edges):
        return None  # impossible
    elif D == 0:
        return []
    for _ in range(10*D):
        candidate_edges = set(random.sample(edges, D))
        # Count valence
        bond_counts = {node: 0 for node in G.nodes()}
        for u, v in G.edges():
            bond_order = 2 if (u, v) in candidate_edges or (v, u) in candidate_edges else 1
            bond_counts[u] += bond_order
            bond_counts[v] += bond_order
        if all(val <= max_valence for val in bond_counts.values()):
            # print("Double bond success")
            return candidate_edges
        else:
            # print('rejecting double bond configuration')
            continue


def graph_to_smiles_with_double(G, double_edges):
    mol = RWMol()
    for _ in G.nodes():
        mol.AddAtom(Chem.Atom("C"))
    for u, v in G.edges():
        btype = Chem.BondType.DOUBLE if (u, v) in double_edges or (v, u) in double_edges else Chem.BondType.SINGLE
        mol.AddBond(int(u), int(v), btype)
    m = mol.GetMol()
    Chem.SanitizeMol(m)
    return Chem.MolToSmiles(m)


def sample_acyclic_hydrocarbons(B, D):
    """Sample n_samples of acyclic hydrocarbons with B total bonds, D of which are double bonds."""
    assert B >= 1, "At least one bond required"
    assert D <= B, "Double bonds must be ≤ total bonds"
    n_carbons = B + 1  # from tree: nodes = edges + 1

    smiles_set = set()
    while len(smiles_set) < 1:
        G = random_tree_prufer(n_carbons, max_degree=4)
        double_bonds = assign_double_bonds(G, D, max_valence=4)
        if double_bonds is None:
            continue
        try:
            smi = graph_to_smiles_with_double(G, double_bonds)
            smiles_set.add(smi)
        except Exception:
            continue  # skip invalid molecules
    return list(smiles_set)[0]

def get_pah_metadata():
    pah_data = pd.read_csv("../data/polycyclic_hydrocarbons/selected_compas_3x.csv")
    # Take the parameters we need for the control
    return pah_data[["bonds", "double_bonds", "atoms"]]

def sample_acyclic_smiles(bond_metadata):
    smiles = bond_metadata.apply(lambda row: 
                                    sample_acyclic_hydrocarbons(row["bonds"], 
                                                                row["double_bonds"]), 
                                                                axis=1)
    return smiles

def write_mols_to_file(ach_df):

    Path("../data/acyclic_hydrocarbons/").mkdir(parents=True, exist_ok=True)

    ach_df["mol"] = ach_df["smiles"].apply(Chem.MolFromSmiles)
    ach_df["filename"] = ach_df.index.to_series().map(lambda x: f"../data/acyclic_hydrocarbons/ach_{x}.mol")
    ach_df.apply(lambda row: Chem.MolToMolFile(row["mol"], row["filename"]), axis=1)
    metadata_file = "../data/acyclic_hydrocarbons/sampled_ach.csv"
    ach_df[["filename", "double_bonds", "bonds", "atoms"]].to_csv(metadata_file)

if __name__ == "__main__":
    # Set random seed for sampling
    random.seed(137)
    # Get controls from ring data
    ach_data = get_pah_metadata()
    # Sample acyclic smiles
    ach_data["smiles"] = sample_acyclic_smiles(ach_data)
    # Write data and meta data to file
    write_mols_to_file(ach_data)
    