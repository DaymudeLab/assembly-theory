import rdkit
from rdkit import Chem

import _pyat
from . import timer

from typing import Optional, Set, Dict, Any
from rdkit import Chem


def molecular_assembly(
    mol: Chem.Mol,
    bounds: Optional[Set[str]] = None,
    no_bounds: bool = False,
    timeout: Optional[int] = None,
    serial: bool = False,
) -> int:
    """
    Computes the molecular assembly index (MA) for a given RDKit molecule.

    Parameters:
    - mol (Chem.Mol): The RDKit molecule object to analyze.
    - bounds (Optional[Set[str]], default=None): A set of bounds to consider. If None, defaults are assigned based on `no_bounds`.
    - no_bounds (bool, default=False): If True, no bounds are used. If False and `bounds` is None, default bounds are applied.
    - timeout (Optional[int], default=None): The time limit in seconds for computation. If exceeded, an error is raised.

    Returns:
    - int: The computed molecular assembly index.

    Raises:
    - ValueError: If `bounds` is specified while `no_bounds` is True.
    - TimeoutError: If computation exceeds the given `timeout`.
    """
    # Convert the molecule to MolBlock format (a string representation).
    mol_block: str = Chem.MolToMolBlock(mol)

    # Validate and initialize bounds.
    bounds = _validate_bounds(bounds, no_bounds)

    if timeout is None:
        # Compute the molecular assembly index with the given bounds.
        ma = _pyat._molecular_assembly(mol_block, bounds, serial)
    else:
        # Run the computation with a timeout to prevent excessive execution time.
        ma = timer.run_with_timeout(
            _pyat._molecular_assembly, timeout, mol_block, bounds, serial
        )

    return ma


def molecular_assembly_verbose(
    mol: Chem.Mol,
    bounds: Optional[Set[str]] = None,
    no_bounds: bool = False,
    timeout: Optional[int] = None,
    serial: bool = False,
) -> Dict[str, int]:
    """
    Computes a verbose molecular assembly index (MA) for a given RDKit molecule,
    returning additional details about the computation.

    Parameters:
    - mol (Chem.Mol): The RDKit molecule object to analyze.
    - bounds (Optional[Set[str]], default=None): A set of bounds to consider. If None, defaults are assigned based on `no_bounds`.
    - no_bounds (bool, default=False): If True, no bounds are used. If False and `bounds` is None, default bounds are applied.
    - timeout (Optional[int], default=None): The time limit in seconds for computation. If exceeded, an error is raised.

    Returns:
    - Dict[str, int]: A dictionary containing:
        - "ma": The computed molecular assembly index.
        - "duplicated_isomorphic_subgraphs": The number of duplicated isomorphic subgraphs.
        - "search_space": The total search space considered.

    Raises:
    - ValueError: If `bounds` is specified while `no_bounds` is True.
    - TimeoutError: If computation exceeds the given `timeout`.
    """
    # Convert the molecule to MolBlock format (a string representation).
    mol_block: str = Chem.MolToMolBlock(mol)

    # Validate and initialize bounds.
    bounds = _validate_bounds(bounds, no_bounds)

    if timeout is None:
        # Compute the verbose molecular assembly index with additional details.
        data = _pyat._molecular_assembly_verbose(mol_block, bounds, serial)
    else:
        # Run the computation with a timeout to prevent excessive execution time.
        data = timer.run_with_timeout(
            _pyat._molecular_assembly_verbose, timeout, mol_block, bounds, serial
        )

    return data


def molecule_info(mol: Chem.Mol) -> str:
    """
    Retrieve molecular information for a given RDKit molecule.

    Parameters:
    - mol (Chem.Mol): An RDKit molecule object.

    Returns:
    - str: a string describing the loaded molecule
    """
    mol_block: str = Chem.MolToMolBlock(mol)  # Convert molecule to MolBlock format.

    info = _pyat._molecule_info(mol_block)  # Extract molecular information.

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
