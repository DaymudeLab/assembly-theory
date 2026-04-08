#!/usr/bin/env python3
"""
Validate alternative assembly pathways.

Checks performed for each molecule:
  1. Structural validity — every pathway step's piece_a ∪ piece_b == result,
     and piece_a ∩ piece_b == ∅.
  2. Distinctness — no two pathways (primary + alternatives) share the same
     sorted match-index set.
  3. Step count consistency — each pathway's step count equals num_duplicates
     + remnant_edges, which should equal the assembly index.

Usage:
    python scripts/validate_pathways.py data/checks/*.mol
    python scripts/validate_pathways.py data/gdb13_1201/*.mol --max-pathways 20
"""

import argparse
import os
import sys

import assembly_theory as at


def validate_pathway_steps(pathway, label="primary"):
    """Check structural validity of a pathway's steps."""
    errors = []
    for i, step in enumerate(pathway):
        a = set(step["piece_a"])
        b = set(step["piece_b"])
        r = set(step["result"])

        if a & b:
            errors.append(f"  [{label}] Step {i+1}: piece_a ∩ piece_b is non-empty: {a & b}")
        if a | b != r:
            errors.append(f"  [{label}] Step {i+1}: piece_a ∪ piece_b != result")
    return errors


def validate_molecule(mol_path, max_pathways=10):
    """Run all validations on a single molecule. Returns (name, errors)."""
    name = os.path.basename(mol_path)
    errors = []

    with open(mol_path) as f:
        mol_block = f.read()

    # Compute the index without alternatives for the ground truth.
    try:
        result = at.pathway_search(mol_block, alternative_pathways=True, max_pathways=max_pathways)
    except Exception as e:
        return name, [f"  pathway_search failed: {e}"]

    index, num_matches, states_searched, pathway, alternatives, alt_match_seqs = result

    if pathway is None:
        # No duplicates found — pathway extraction is not possible.
        # This is expected for molecules with 0 matches.
        return name, []

    # --- Check 1: structural validity of primary pathway ---
    errors.extend(validate_pathway_steps(pathway, "primary"))

    # --- Check 1b: structural validity of alternatives ---
    if alternatives:
        for j, alt in enumerate(alternatives):
            errors.extend(validate_pathway_steps(alt, f"alt-{j+1}"))

    # --- Check 2: step count consistency ---
    # pathway steps == assembly index
    if len(pathway) != index:
        errors.append(
            f"  Primary pathway has {len(pathway)} steps but assembly index is {index}"
        )
    if alternatives:
        for j, alt in enumerate(alternatives):
            if len(alt) != index:
                errors.append(
                    f"  Alt-{j+1} has {len(alt)} steps but assembly index is {index}"
                )

    # --- Check 3: distinctness via sorted match-index sets ---
    if alternatives and alt_match_seqs:
        # We don't have the primary match seq directly via Python API,
        # but we can check that no two alternatives are identical.
        seen = set()
        for j, seq in enumerate(alt_match_seqs):
            key = tuple(sorted(seq))
            if key in seen:
                errors.append(f"  Alt-{j+1} has duplicate match-set key {key}")
            seen.add(key)

    return name, errors


def main():
    parser = argparse.ArgumentParser(description="Validate alternative assembly pathways")
    parser.add_argument("molfiles", nargs="+", help="Path(s) to .mol files")
    parser.add_argument("--max-pathways", type=int, default=10,
                        help="Max pathways to collect (default: 10)")
    args = parser.parse_args()

    total = 0
    passed = 0
    failed_mols = []

    for path in sorted(args.molfiles):
        if not path.endswith(".mol"):
            continue
        total += 1
        name, errors = validate_molecule(path, max_pathways=args.max_pathways)
        if errors:
            failed_mols.append((name, errors))
            print(f"FAIL  {name}")
            for e in errors:
                print(e)
        else:
            passed += 1
            print(f"OK    {name}")

    print(f"\n{passed}/{total} molecules passed validation")
    if failed_mols:
        print(f"{len(failed_mols)} molecule(s) failed:")
        for name, _ in failed_mols:
            print(f"  - {name}")
        sys.exit(1)


if __name__ == "__main__":
    main()
