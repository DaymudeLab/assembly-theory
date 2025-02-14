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
    
    bounds = _validate_bounds(bounds, no_bounds)
    
    ma = _compute_index(mol_block, bounds)  # Compute the molecular index using the given bounds.
    
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

    bounds = _validate_bounds(bounds, no_bounds)

    data = _compute_verbose_index(mol_block, bounds)  # Compute verbose molecular index.

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

    info = _molecule_info(mol_block)  # Extract molecular information.
    
    return info

def _validate_bounds(bounds: Optional[Set[str]], no_bounds: bool) -> Set[str]:
    """
    Validates and initializes the `bounds` variable based on `no_bounds` flag.

    Args:
        bounds (Optional[Set[str]]): The initial bounds, if any.
        no_bounds (bool): Flag indicating whether bounds should be absent.

    Returns:
        Set[str]: A set containing bounds if applicable.

    Raises:
        ValueError: If `bounds` is specified but `no_bounds` is True.
    """
    if bounds is None:
        if no_bounds:
            return set()  # Initialize an empty set if no bounds are provided.
        else:
            return {"addition", "vector"}
    elif (bounds is not None) and no_bounds:
        raise ValueError("bounds specified but `no_bounds` is True.")
    
    return bounds
