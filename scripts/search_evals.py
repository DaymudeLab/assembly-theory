import assembly_theory as at
import matplotlib.pyplot as plt
import os
import argparse
import re

# CLI parsing
parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Path to directory to be benchmarked")
parser.add_argument("--progress", action='store_true', help="Print progress of benchmark")
parser.add_argument("--out", help="Directory to place output files")

args = parser.parse_args()

dir = args.directory
if dir[-1] != '/':
    dir += '/'
num_dups = []
states_searched = []

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
        mol = mol_file.read()
        (_, dup, state) = at.index_search(
            mol,
            enumerate_str = 'grow-erode',
            canonize_str = 'tree-nauty',
            parallel_str = 'depth-one',
            memoize_str = 'canon-index',
            kernel_str = 'none',
            bound_strs = {'int', 'vec-simple', 'vec-small-frags'}
        )
        if dup > 2500:
            continue
        num_dups.append(dup)
        states_searched.append(state)
    
    if args.progress:
        print(str(i) + "/" + str(total))

# Plotting
plt.scatter(x=num_dups, y=states_searched)
plt.xlabel("Num duplicatable matches")
plt.ylabel("States Searched")

# Save eval data and figure
regex = re.compile('/?([^/]*)/$')
data_dir = regex.search(dir).group(1)
out_dir = args.out if args.out else "./"
if out_dir[-1] != "/":
    out_dir += "/"

os.makedirs(out_dir, exist_ok=True)
plt.savefig(out_dir + data_dir + ".svg", format="svg")

