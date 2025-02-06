import rdkit
from rdkit import Chem
from .orca import assembly_index
def compute_ma(mol):

    mol_block = rdkit.Chem.MolToMolBlock(mol)
    ma = assembly_index(mol_block)

    return ma
