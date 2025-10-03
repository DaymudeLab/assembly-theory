import requests
import pandas as pd
from pathlib import Path

import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def download_compas_csv(save_path="compas-3x.csv"):
    url = "https://gitlab.com/porannegroup/compas/-/raw/32b1a7904ebc525ec0095edc41daba0b9a2b3436/COMPAS-3/compas-3x.csv"
    response = requests.get(url)
    
    if response.status_code == 200:
        with open(save_path, "wb") as f:
            f.write(response.content)
        print(f"Downloaded to {save_path}")
    else:
        raise Exception(f"Failed to download file: {response.status_code}")

def count_double_bonds(mol):
    rdkit.Chem.Kekulize(mol)
    double_bond_count = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            double_bond_count += 1
    return double_bond_count

def read_and_reduce_compas_3x():
    compas_df = pd.read_csv("compas-3x.csv")
    compas_df = compas_df.query("n_rings <= 8")
    compas_df["mol"] = compas_df["smiles"].apply(Chem.MolFromSmiles)
    compas_df["formula"] = compas_df["mol"].apply(Chem.rdMolDescriptors.CalcMolFormula)
    compas_df["bonds"] = [c.GetNumBonds() for c in compas_df["mol"]]
    compas_df["double_bonds"] = compas_df["mol"].apply(count_double_bonds)
    compas_df["atoms"] = [c.GetNumAtoms() for c in compas_df["mol"]]

    reduced_df = compas_df[["smiles",
                        "mol",
                        "n_rings",
                        "double_bonds",
                        "bonds",
                        "atoms"]].reset_index()
    return reduced_df

def write_mols_to_file(compas_df):
    Path("../data/polycyclic_hydrocarbons/").mkdir(parents=True, exist_ok=True)

    compas_df["filename"] = compas_df.index.to_series().map(lambda x: f"../data/polycyclic_hydrocarbons/pah_{x}.mol")
    compas_df.apply(lambda row: Chem.MolToMolFile(row["mol"], row["filename"]), axis=1)
    meta_data_file = "../data/polycyclic_hydrocarbons/selected_compas_3x.csv"
    compas_df[["filename", "n_rings", "double_bonds", "bonds", "atoms"]].to_csv(meta_data_file)

if __name__ == "__main__":
    # Grab COMPAS-3 Data from gitlab 
    # Cite as A. Wahab and R. Gershoni-Poranne, COMPAS-3, DOI: 10.1039/D4CP01027B
    download_compas_csv()
    
    # Read and reduce data
    reduced_compas = read_and_reduce_compas_3x()

    # Generate filenames for mol and write to file
    write_mols_to_file(reduced_compas)