import assembly_theory as at
import os
import argparse
import re
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

# Function to run benchmark
def run_bench(dir):
    df = pd.DataFrame(columns=["File", "Atoms", "Bonds", "Rings", "Index", "Dupes", "StatesSearched"])

    if dir[-1] != '/':
        dir += '/'

    # Find number of files to be benchmarked
    total = 0
    if args.progress:
        for file in os.listdir(os.fsencode(dir)):
            filename = os.fsdecode(file)
            
            if filename.endswith('.mol'):
                total += 1

    # Loop through mol files and compute MA
    for (i, file) in enumerate(os.listdir(os.fsencode(dir))):
        filename = os.fsdecode(file)
        
        if not filename.endswith('.mol'):
            continue

        with open(dir + filename) as mol_file:
            mol_block = mol_file.read()
            mol = Chem.MolFromMolBlock(mol_block)

            (ma, dup, state) = at.index_search(
                mol_block,
                enumerate_str = 'grow-erode',
                canonize_str = 'tree-nauty',
                parallel_str = 'none',
                memoize_str = 'canon-index',
                kernel_str = 'none',
                bound_strs = ['int', 'vec-simple', 'vec-small-frags']
            )
            
            new_row = pd.DataFrame({
                'File': [filename],
                'Atoms': [mol.GetNumAtoms()],
                'Bonds': [mol.GetNumBonds()],
                'Rings': [rdMolDescriptors.CalcNumRings(mol)],
                'Index': [ma],
                'Dupes': [dup],
                'StatesSearched': [state],
            })
            df = pd.concat([df, new_row])
        
        if args.progress:
            print(str(i) + "/" + str(total))

    return df


# CLI parsing
parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Path to directory to be benchmarked or file to be read")
parser.add_argument("--progress", action='store_true', help="Print progress of benchmark")
parser.add_argument("--out", help="Directory to place output files")
parser.add_argument("--name", help="Evaluation name")

# Parse arguments
args = parser.parse_args()
dir = args.directory

out_dir = args.out if args.out else "./"
if out_dir[-1] != "/":
    out_dir += "/"
os.makedirs(out_dir, exist_ok=True)

data_name = ""
if args.name:
    # Data name is provided by user
    data_name = args.name
else:
    # Get data name from benchmark directory name
    regex = re.compile('/?([^/]*)/$')
    data_name = regex.search(dir).group(1)

out_file = out_dir + data_name + '.out'

# Fill lists of molecule assemlby info
df = run_bench(dir)
df.to_csv(out_file, index=False)
