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
    raise Exception("Invald input type")

out_dir = args.out if args.out else os.path.dirname(file_str)
if out_dir[-1] != "/":
    out_dir += "/"
os.makedirs(out_dir, exist_ok=True)
data_name = args.name if args.name else os.path.basename(file_str)[:-4]
out_file = out_dir + data_name + '.svg'

# Plotting
# dataframe columns: File, Atoms, Bonds, Rings, Index, Dupes, StatesSearched
df = pd.read_csv(args.file)

plt.scatter(x=df['Dupes'], y=df['StatesSearched'], c=df['Bonds'])
plt.xlabel('Num duplicatable matches')
plt.ylabel('States Searched')
plt.savefig(out_file, format='svg')