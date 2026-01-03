#!/usr/bin/env python3
"""
Validate and copy a custom species tree for use in TimeTree pipeline.

This script validates that:
1. The tree file exists and is readable
2. The tree is in valid Newick format
3. The tree is properly rooted (bifurcating root)
4. The leaf names exactly match species in the supermatrix
"""

import sys
from pathlib import Path
from Bio import SeqIO
from ete3 import Tree


def get_species_from_fasta(fasta_file: Path) -> set:
    """Extract species names from FASTA file."""
    species = set()
    for record in SeqIO.parse(fasta_file, "fasta"):
        species.add(record.id)
    return species


def validate_tree_rooted(tree: Tree) -> bool:
    """
    Check if tree is properly rooted (bifurcating root).

    A properly rooted tree has exactly 2 children at the root.
    """
    root_children = tree.get_children()
    return len(root_children) == 2


def validate_species_match(tree_species: set, expected_species: set) -> tuple:
    """
    Validate that tree species match expected species.

    Returns:
        (is_valid, missing_from_tree, extra_in_tree)
    """
    missing_from_tree = expected_species - tree_species
    extra_in_tree = tree_species - expected_species
    is_valid = len(missing_from_tree) == 0 and len(extra_in_tree) == 0
    return is_valid, missing_from_tree, extra_in_tree


def main():
    # Get input/output paths from Snakemake
    tree_file = Path(snakemake.input.tree)
    supermatrix_file = Path(snakemake.input.supermatrix)
    out_file = Path(snakemake.output.rooted)
    log_file = Path(snakemake.log[0])

    # Setup logging
    log_file.parent.mkdir(parents=True, exist_ok=True)

    errors = []

    with open(log_file, 'w') as log:
        log.write(f"Validating custom species tree: {tree_file}\n")
        log.write("=" * 60 + "\n\n")

        # Step 1: Check tree file exists
        if not tree_file.exists():
            msg = f"ERROR: Tree file not found: {tree_file}"
            log.write(msg + "\n")
            errors.append(msg)
        else:
            log.write(f"[OK] Tree file exists: {tree_file}\n")

        # Step 2: Try to parse tree
        try:
            tree = Tree(str(tree_file))
            log.write(f"[OK] Tree parsed successfully (Newick format valid)\n")
        except Exception as e:
            msg = f"ERROR: Failed to parse tree file as Newick: {e}"
            log.write(msg + "\n")
            errors.append(msg)
            tree = None

        if tree:
            # Step 3: Check if tree is rooted
            if validate_tree_rooted(tree):
                log.write(f"[OK] Tree is properly rooted (bifurcating root)\n")
            else:
                root_children = len(tree.get_children())
                msg = (
                    f"ERROR: Tree is not properly rooted.\n"
                    f"       Root has {root_children} children (expected 2 for bifurcating root).\n"
                    f"       Please provide a rooted tree or use tree_method='concatenation' or 'coalescent'."
                )
                log.write(msg + "\n")
                errors.append(msg)

            # Step 4: Validate species names
            tree_species = set(tree.get_leaf_names())
            expected_species = get_species_from_fasta(supermatrix_file)

            log.write(f"\nSpecies in tree: {len(tree_species)}\n")
            log.write(f"Species in data: {len(expected_species)}\n\n")

            is_valid, missing, extra = validate_species_match(tree_species, expected_species)

            if is_valid:
                log.write(f"[OK] All species names match between tree and data\n")
            else:
                if missing:
                    msg = f"ERROR: Species missing from tree (present in data): {sorted(missing)}"
                    log.write(msg + "\n")
                    errors.append(msg)
                if extra:
                    msg = f"ERROR: Extra species in tree (not in data): {sorted(extra)}"
                    log.write(msg + "\n")
                    errors.append(msg)

        # Final status
        log.write("\n" + "=" * 60 + "\n")
        if errors:
            log.write(f"FAILED: {len(errors)} validation error(s)\n")
            for i, err in enumerate(errors, 1):
                log.write(f"  {i}. {err}\n")
            raise ValueError(
                f"Custom tree validation failed with {len(errors)} error(s).\n"
                f"See log file for details: {log_file}\n"
                + "\n".join(errors)
            )
        else:
            log.write("SUCCESS: All validations passed\n")

            # Write validated tree to output
            out_file.parent.mkdir(parents=True, exist_ok=True)
            tree.write(outfile=str(out_file), format=1)
            log.write(f"\nTree copied to: {out_file}\n")

            print(f"Custom tree validated successfully: {len(tree_species)} species")
            print(f"Output: {out_file}")


if __name__ == "__main__":
    main()
