import rdkit
from rdkit import Chem

from .orca import *
from .orca import __all__
from .orca import __doc__

def compute_ma(mol, bounds=None, no_bounds=False):

    mol_block = rdkit.Chem.MolToMolBlock(mol)
    ma = compute_index(mol_block)

    return ma

def compute_ma_verbose(mol, bounds=None, no_bounds=False):
    mol_block = rdkit.Chem.MolToMolBlock(mol)
    data = compute_verbose_index(mol_block)

    return data

def get_molecule_info(mol):
    mol_block = rdkit.Chem.MolToMolBlock(mol)
    info = molecule_info(mol_block)
    return info
