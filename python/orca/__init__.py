import rdkit
from rdkit import Chem

from .orca import *
from .orca import __all__
from .orca import __doc__

from typing import Optional, Set, Dict, Any
from rdkit import Chem

def compute_ma(mol: Chem.Mol, bounds: Optional[Set[str]] = None, no_bounds: bool = False) -> int:
    """
    Compute the molecular index (MA) for a given RDKit molecule.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.
    - bounds (Optional[Set[str]]): A set of bounds to consider (default: None).
    - no_bounds (bool): Flag to indicate whether bounds should be ignored (default: False).

    Returns:
    - int: The computed assembly index.
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format (string representation).
    
    if (bounds is None) and no_bounds:
        bounds = set()  # Initialize an empty set if no bounds are provided.
    elif (bounds is None) and not no_bounds:
        bounds = set(["addition", "vector"])
    elif (bounds is not None) and not no_bounds:
        raise ValueError("bounds specified but `no_bound` is True.")
    
    ma = compute_index(mol_block, bounds)  # Compute the molecular index using the given bounds.
    
    return ma


def compute_ma_verbose(mol: Chem.Mol, bounds: Optional[Set[str]] = None, no_bounds: bool = False) -> Dict[str, int]:
    """
    Compute the verbose molecular index, returning additional details.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.
    - bounds (Optional[Set[str]]): A set of bounds to consider (default: None).
    - no_bounds (bool): Flag to indicate whether bounds should be ignored (default: False).

    Returns:
    - Dict[str, int]: A dictionary containing the assembly index, the number of duplicated isomorphic subgraphs, and the search space.
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format (string representation).

    if (bounds is None) and no_bounds:
        bounds = set()  # Initialize an empty set if no bounds are provided.
    elif (bounds is None) and not no_bounds:
        bounds = set(["addition", "vector"])
    elif (bounds is not None) and not no_bounds:
        raise ValueError("bounds specified but `no_bound` is True.")
    data = compute_verbose_index(mol_block, bounds)  # Compute verbose molecular index.

    return data


def get_molecule_info(mol: Chem.Mol) -> str:
    """
    Retrieve molecular information for a given RDKit molecule.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.

    Returns:
    - str: a string describing the loaded molecule
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format.

    info = molecule_info(mol_block)  # Extract molecular information.
    
    return info

