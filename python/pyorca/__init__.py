import rdkit
from rdkit import Chem

from . import _pyorca
from ._pyorca import __all__
from ._pyorca import __doc__
from . import timer

from typing import Optional, Set, Dict, Any
from rdkit import Chem

def compute_ma(mol: Chem.Mol, 
               bounds: Optional[Set[str]] = None, 
               no_bounds: bool = False,
               timeout: Optional[int] = None) -> int:
    """
    Compute the molecular index (MA) for a given RDKit molecule.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.
    - bounds (Optional[Set[str]]): A set of bounds to consider (default: None).
    - no_bounds (bool): Flag to indicate whether bounds should be ignored (default: False).
    - timeout: (Optional[int]): time in seconds to allow computation to proceed, if timeout is exceed function will be terminated and error will be raised.
    
    Returns:
    - int: The computed assembly index.
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format (string representation).
    
    bounds = _validate_bounds(bounds, no_bounds)
    
    if timeout is None:
        ma = _compute_index(mol_block, bounds)  # Compute the molecular index using the given bounds.
    else:
        ma = timer.run_with_timeout(_compute_index, timeout, mol_block, bounds)
    
    return ma


def compute_ma_verbose(mol: Chem.Mol, 
                       bounds: Optional[Set[str]] = None, 
                       no_bounds: bool = False, 
                       timeout: Optional[int] = None) -> Dict[str, int]:
    """
    Compute the verbose molecular index, returning additional details.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.
    - bounds (Optional[Set[str]]): A set of bounds to consider (default: None).
    - no_bounds (bool): Flag to indicate whether bounds should be ignored (default: False).
    - timeout: (Optional[int]): time in seconds to allow computation to proceed, if timeout is exceed function will be terminated and error will be raised.
    
    Returns:
    - Dict[str, int]: A dictionary containing the assembly index, the number of duplicated isomorphic subgraphs, and the search space.
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format (string representation).

    bounds = _validate_bounds(bounds, no_bounds)
    
    if timeout is None:
        data = _compute_verbose_index(mol_block, bounds)  # Compute the molecular index using the given bounds.
    else:
        data = timer.run_with_timeout(_compute_verbose_index, timeout, mol_block, bounds)
    
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
            return {"intchain", "vecchain"}
    elif (bounds is not None) and no_bounds:
        raise ValueError("bounds specified but `no_bounds` is True.")
    
    return bounds
