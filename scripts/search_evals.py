import assembly_theory as at
import matplotlib.pyplot as plt
import os
import argparse
import re

# Function to run benchmark
def run_bench(dir, indices, num_dups, states_searched):
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
            mol = mol_file.read()
            (ma, dup, state) = at.index_search(
                mol,
                enumerate_str = 'grow-erode',
                canonize_str = 'tree-nauty',
                parallel_str = 'depth-one',
                memoize_str = 'canon-index',
                kernel_str = 'none',
                bound_strs = {'int', 'vec-simple', 'vec-small-frags'}
            )
            indices.append(ma)
            num_dups.append(dup)
            states_searched.append(state)
        
        if args.progress:
            print(str(i) + "/" + str(total))


# Function for reading from file
def read_from_file(file, indices, num_dups, states_searched):
    with open(file) as f:
        f.readline()

        # Read assembly indices
        line = f.readline()
        while line != "Num duplicate matches\n":
            indices.append(int(line.strip()))
            line = f.readline()
            
        
        # Read num duplicates
        line = f.readline()
        while line != "States Searched\n":
            num_dups.append(int(line.strip()))
            line = f.readline()
        
        # Read states searched
        line = f.readline()
        while line:
            states_searched.append(int(line.strip()))
            line = f.readline()


# CLI parsing
parser = argparse.ArgumentParser()
parser.add_argument("directory", help="Path to directory to be benchmarked or file to be read")
parser.add_argument("--progress", action='store_true', help="Print progress of benchmark")
parser.add_argument("--out", help="Directory to place output files")
parser.add_argument("--name", help="Evaluation name")

# Parse arguments
args = parser.parse_args()
dir = args.directory

reading = dir[-4:] == ".out"

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
    if reading:
        regex = re.compile("/([^/]*)\\.out$")
    else:
        regex = re.compile('/?([^/]*)/$')
    data_name = regex.search(dir).group(1)

out_file = out_dir + data_name

# Fill lists of molecule assemlby info
indices = []
num_dups = []
states_searched = []

if reading:
    read_from_file(dir, indices, num_dups, states_searched)
else:
    run_bench(dir, indices, num_dups, states_searched)

# Plotting
plt.scatter(x=num_dups, y=states_searched)
plt.xlabel("Num duplicatable matches")
plt.ylabel("States Searched")

# Save eval figure
plt.savefig(out_file + ".svg", format="svg")

# Save eval data if benchmarking
if not reading:
    with open(out_file + ".out", 'w') as f:
        f.write("Assembly indices\n")
        for x in indices:
            f.write(str(x) + '\n')

        f.write("Num duplicate matches\n")
        for x in num_dups:
            f.write(str(x) + '\n')
        
        f.write("States Searched\n")
        for x in states_searched:
            f.write(str(x) + '\n')
