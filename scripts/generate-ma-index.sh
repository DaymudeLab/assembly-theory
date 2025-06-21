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
select exec_choice in "assembly_go (Jirasek et al., 2024)" "assembly_cpp (Seet et al., 2024)" "assembly-theory"
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
            if [ ! -f "assembly_cpp" ]; then
                echo -n "ERROR: Missing ./assembly_cpp executable "
                echo "(provided privately by Seet et al.)."
                exit 1
            fi
            executable="./assembly_cpp"
            break
            ;;
        3)
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

    # assembly_go and assembly-theory expect "<molecule>.mol" but assembly_cpp
    # expects only "<molecule>" with the ".mol" part stripped off. Also,
    # assembly_cpp prints a ton of unnecessary information requiring some
    # parsing to get just the assembly index. Lastly, assembly_cpp generates
    # auxiliary output files that need to be removed.
    if [ $executable = "./assembly_cpp" ]; then
        molpath_stripped=$(echo "$direntry" | sed -e "s/.mol//g")
        maindex=$("$executable" "$molpath_stripped" -pathway=0 | tail -n 1 | awk '{print $NF}')
        rm "${molpath_stripped}Out"
    else
        maindex=$("$executable" "$direntry")
    fi
    echo "$molfile,$maindex" >> "$mafile"
done
echo ""
