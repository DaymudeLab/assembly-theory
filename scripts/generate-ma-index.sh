#!/bin/bash

# Get target dataset.
PS3="Generate ma-index.csv for: "
select dataset in ../data/*
do
    case $dataset in
        "")
            echo "ERROR: Invalid option $REPLY."
            ;;
        *)
            break
            ;;
    esac
done

# Let the user choose which executable should generate ground truth.
PS3="Calculate assembly indices using: "
select exec_choice in "assembly_go" "assembly-theory"
do
    case $REPLY in
        1)
            if [ ! -f "assembly_go" ]; then
                echo -n "ERROR: Missing ./assembly_go executable "
                echo "(https://github.com/croningp/assembly_go)."
                exit 1
            fi
            executable="./assembly_go"
            break
            ;;
        2)
            echo "Building a release version of assembly-theory..."
            cargo build --release
            executable="./../target/release/assembly-theory"
            break
            ;;
        *)
            echo "ERROR: Invalid option $REPLY."
            ;;
    esac
done

# Initialize the ma-index.csv file.
mafile="$dataset/ma-index.csv"
> "$mafile"
echo "file_name,assembly_idx" >> "$mafile"

# Calculate and record assembly index for all .mol files in the dataset.
for direntry in "$dataset"/*.mol
do
    molfile=$(basename "$direntry")
    echo -ne "\r\e[K$exec_choice: Calculating assembly index of $molfile..."
    maindex=$("$executable" "$direntry")
    echo "$molfile,$maindex" >> "$mafile"
done
echo ""
