import assembly_theory as at
import rdkit.Chem
import pytest 

def test_timeout_pass():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly_verbose(mol, timeout=5.0) == {'duplicates': 20, 'index': 8, 'space': 36}

    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert at.molecular_assembly(mol, timeout=5.0) == 8

def test_timeout_fail():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    with pytest.raises(TimeoutError):
       at.molecular_assembly_verbose(mol, timeout=0.0)

    with pytest.raises(TimeoutError):
       at.molecular_assembly(mol, timeout=0.0)
