import pyorca
import rdkit.Chem
import pytest 

def test_timeout_pass():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert pyorca.compute_ma_verbose(mol, timeout=1.0) == {'duplicates': 20, 'index': 8, 'space': 36}

    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    assert pyorca.compute_ma(mol, timeout=1.0) == 8

def test_timeout_fail():
    aspirin_smi = "O=C(C)Oc1ccccc1C(=O)O"
    mol = rdkit.Chem.MolFromSmiles(aspirin_smi)
    with pytest.raises(TimeoutError):
       pyorca.compute_ma_verbose(mol, timeout=0.0)

    with pytest.raises(TimeoutError):
       pyorca.compute_ma(mol, timeout=0.0)
