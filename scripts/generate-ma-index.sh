#!/bin/bash

# Get target dataset and check that it exists.
read -p "Generate ma-index.csv for: assembly-theory/data/" dataset
if [ ! -d "../data/$dataset" ]; then
    echo "ERROR: assembly-theory/data/$dataset does not exist."
fi
datadir="../data/$dataset"

# Build a release version of assembly-theory.
echo "Building a release version of assembly-theory..."
cargo build --release

# Initialize the ma-index.csv file.
mafile="$datadir/ma-index.csv"
> "$mafile"
echo "file_name,assembly_idx" >> "$mafile"

# Calculate and record assembly index for all .mol files in the dataset.
for direntry in "$datadir"/*.mol
do
    molfile=$(basename "$direntry")
    echo -ne "\r\e[KCalculating assembly index of $molfile..."
    maindex=$(../target/release/assembly-theory "$direntry")
    echo "$molfile,$maindex" >> "$mafile"
done
echo ""
