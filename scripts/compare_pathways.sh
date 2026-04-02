#!/usr/bin/env zsh
# =============================================================================
# Compare assembly pathway outputs between Rust and C++ implementations.
#
# BACKGROUND
# ----------
# Both codebases compute a molecule's assembly index by finding pairs of
# edge-disjoint isomorphic subgraphs ("duplicates") in a top-down search.
# The assembly index = total_edges - 1 - sum(dup_sizes - 1).
#
# OUTPUT FORMATS
# --------------
# Rust (--extract-pathway):
#   Prints a bottom-up pathway (explicit join steps) and a decomposition
#   summary with duplicate pairs as vertex-pair edge lists.
#   "Remnant edges" = total_edges - sum(duplicate Right sizes)
#     i.e., all edges not accounted for by duplicate removal.
#
# C++ (writes <name>Pathway JSON):
#   Outputs the molecule graph, a "remnant" subgraph, a "duplicates" list,
#   and a "removed_edges" list.
#   The C++ preprocesses the molecule BEFORE search: edges whose chemical
#   signature (atom_type + bond_type + atom_type) appears only once in the
#   molecule are stripped as "removed_edges" — they can never be part of a
#   duplicate pair. The C++ "remnant" only contains the remaining edges.
#   So: C++ remnant + len(removed_edges) = Rust remnant.
#
# WHAT THIS SCRIPT COMPARES
# -------------------------
# For each .mol file in data/checks/, runs both implementations and compares:
#   1. Assembly index (should be identical)
#   2. Number of duplicate pairs (should be identical)
#   3. Duplicate sizes, sorted descending (should be identical)
#   4. Remnant edge count: Rust remnant vs C++ (remnant + removed_edges)
#      These should match because both equal total_edges - sum(dup sizes).
#
# The specific edges in the remnant/duplicates may differ between runs (the
# search is nondeterministic due to parallelism and HashMap ordering), but the
# structural invariants above must always agree.
#
# NOTES
# -----
# - The C++ binary expects input WITHOUT .mol extension; it appends .mol
#   internally. It also requires the .mol file in the current directory.
# - Some molecules (morphine, cyanidin) may take minutes in the C++ binary.
# - Expects assemblycpp-v5 to be checked out alongside this repository
#   (i.e., at ../../assemblycpp-v5 relative to this script).
# =============================================================================

set -euo pipefail

SCRIPT_DIR="${0:A:h}"
RUST_DIR="$SCRIPT_DIR/.."
CPP_DIR="$SCRIPT_DIR/../../assemblycpp-v5"
RUST_BIN="$RUST_DIR/target/debug/assembly-theory"
CPP_BIN="$CPP_DIR/build/bin/assembly"
MOL_DIR="$RUST_DIR/data/checks"
WORK_DIR="$(mktemp -d)"

trap 'rm -rf "$WORK_DIR"' EXIT

# Build Rust binary if needed.
if [[ ! -x "$RUST_BIN" ]]; then
    echo "Building Rust binary..."
    (cd "$RUST_DIR" && cargo build)
fi

# Check C++ binary exists.
if [[ ! -x "$CPP_BIN" ]]; then
    echo "ERROR: C++ binary not found at $CPP_BIN"
    echo "Build it first: cd $CPP_DIR/build && cmake .. && make"
    exit 1
fi

pass=0
fail=0
errors=()

for molfile in "$MOL_DIR"/*.mol; do
    molname="${molfile:t:r}"  # basename without extension
    echo -n "Testing $molname ... "

    # --- Rust ---
    rust_output=$("$RUST_BIN" --extract-pathway "$molfile" 2>&1)
    rust_index=$(echo "$rust_output" | head -1)
    rust_dups=$(echo "$rust_output" | grep "Duplicates:" | sed 's/.*Duplicates: \([0-9]*\) (sizes: \[\(.*\)\])/\1|\2/')
    rust_num_dups=$(echo "$rust_dups" | cut -d'|' -f1)
    rust_dup_sizes=$(echo "$rust_dups" | cut -d'|' -f2)
    rust_remnant=$(echo "$rust_output" | grep "Remnant edges:" | sed 's/.*Remnant edges: //')

    # --- C++ ---
    # C++ needs mol file in cwd, named without extension.
    cp "$molfile" "$WORK_DIR/${molname}.mol"
    (cd "$WORK_DIR" && "$CPP_BIN" "$molname" > /dev/null 2>&1)

    cpp_out_file="$WORK_DIR/${molname}Out"
    cpp_pathway_file="$WORK_DIR/${molname}Pathway"

    if [[ ! -f "$cpp_out_file" ]]; then
        echo "SKIP (C++ produced no output)"
        continue
    fi

    cpp_index=$(grep "assembly index:" "$cpp_out_file" | sed 's/.*assembly index: //')

    # Parse the C++ pathway JSON with python for robustness.
    if [[ ! -f "$cpp_pathway_file" ]]; then
        echo "SKIP (C++ produced no pathway)"
        continue
    fi

    cpp_summary=$(python3 -c "
import json, sys
with open('$cpp_pathway_file') as f:
    data = json.load(f)
dups = data.get('duplicates', [])
sizes = sorted([len(d['Left']) for d in dups], reverse=True)
remnant_edges = len(data.get('remnant', [{}])[0].get('Edges', []))
removed_edges = len(data.get('removed_edges', []))
print(f'{len(dups)}|{\",\".join(str(s) for s in sizes)}|{remnant_edges + removed_edges}')
")
    cpp_num_dups=$(echo "$cpp_summary" | cut -d'|' -f1)
    cpp_dup_sizes=$(echo "$cpp_summary" | cut -d'|' -f2)
    cpp_remnant=$(echo "$cpp_summary" | cut -d'|' -f3)

    # --- Compare ---
    # Sort duplicate sizes for comparison (both should be descending already).
    rust_sizes_sorted=$(echo "$rust_dup_sizes" | tr ',' '\n' | tr -d ' ' | sort -rn | tr '\n' ',' | sed 's/,$//')
    cpp_sizes_sorted=$(echo "$cpp_dup_sizes" | tr ',' '\n' | tr -d ' ' | sort -rn | tr '\n' ',' | sed 's/,$//')

    ok=true
    details=""

    if [[ "$rust_index" != "$cpp_index" ]]; then
        ok=false
        details+="  index: rust=$rust_index cpp=$cpp_index\n"
    fi
    if [[ "$rust_num_dups" != "$cpp_num_dups" ]]; then
        ok=false
        details+="  num_dups: rust=$rust_num_dups cpp=$cpp_num_dups\n"
    fi
    if [[ "$rust_sizes_sorted" != "$cpp_sizes_sorted" ]]; then
        ok=false
        details+="  dup_sizes: rust=[$rust_sizes_sorted] cpp=[$cpp_sizes_sorted]\n"
    fi
    if [[ "$rust_remnant" != "$cpp_remnant" ]]; then
        ok=false
        details+="  remnant: rust=$rust_remnant cpp=$cpp_remnant\n"
    fi

    if $ok; then
        echo "PASS (index=$rust_index dups=$rust_num_dups sizes=[$rust_sizes_sorted] remnant=$rust_remnant)"
        pass=$((pass + 1))
    else
        echo "FAIL"
        echo -e "$details"
        errors+=("$molname")
        fail=$((fail + 1))
    fi
done

echo ""
echo "=== Results: $pass passed, $fail failed ==="
if (( ${#errors[@]} > 0 )); then
    echo "Failed molecules: ${errors[*]}"
    exit 1
fi

# Clean up any C++ output files left from manual runs in the C++ directory.
for molfile in "$MOL_DIR"/*.mol; do
    molname="${molfile:t:r}"
    rm -f "$CPP_DIR/${molname}Out" "$CPP_DIR/${molname}Pathway" \
          "$CPP_DIR/${molname}.mol" "$CPP_DIR/${molname}.molOut" \
          "$CPP_DIR/${molname}.molPathway"
done
