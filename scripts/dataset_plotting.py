import matplotlib.pyplot as plt
import pandas as pd
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument("file", help=".out file produced by dataset_info.py")
parser.add_argument("--out", help="Directory to place output files")
parser.add_argument("--name", help="Evaluation name")

# Parse arguments
args = parser.parse_args()
file_str = args.file
if len(file_str) < 4 or file_str[-4:] != '.out':
    raise Exception("Invaldi input type")

out_dir = args.out if args.out else os.path.dirname(file_str)
data_name = args.name if args.name else os.path.basename(file_str)
out_file = out_dir + data_name + '.svg'

# Plotting
df = pd.read_csv(args.file)

plt.scatter(x=df['Dupes'], y=df['StatesSearched'])
plt.xlabel('Num duplicatable matches')
plt.ylabel('States Searched')
plt.savefig(out_file, format='svg')