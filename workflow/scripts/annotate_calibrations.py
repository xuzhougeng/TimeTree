#!/usr/bin/env python3
"""
Annotate rooted tree with fossil/node calibrations for MCMCtree.
Reads calibration table (TSV) and injects constraints at MRCA nodes.

MCMCtree calibration syntax (placed before node in Newick):
  '>0.5<1.0'  - min 0.5, max 1.0 (uniform/truncated bounds)
  'B(0.5,1.0)' - soft bounds with 2.5% tail probabilities
  'L(0.5)'    - lower bound only
  'U(1.0)'    - upper bound only
"""

import sys
from pathlib import Path
import pandas as pd
from ete3 import Tree


def format_calibration(min_age: float, max_age: float, prior: str = None) -> str:
    """
    Format calibration constraint for MCMCtree.
    Default: B(min,max) soft bounds.
    """
    if prior and prior.upper().startswith('B'):
        return f"'B({min_age},{max_age})'"
    elif prior and prior.upper().startswith('L'):
        return f"'L({min_age})'"
    elif prior and prior.upper().startswith('U'):
        return f"'U({max_age})'"
    else:
        # Default: soft bounds B(min,max)
        return f"'B({min_age},{max_age})'"


def find_mrca_and_annotate(tree: Tree, taxa_list: "list[str]", calibration_str: str, label: str = "") -> "tuple[bool, str]":
    """
    Find MRCA of taxa list and add calibration annotation.
    Returns (success, message).
    """
    leaf_names = set(tree.get_leaf_names())

    # Validate taxa
    valid_taxa = [t for t in taxa_list if t in leaf_names]
    missing = [t for t in taxa_list if t not in leaf_names]

    if missing:
        return False, f"Taxa not found in tree: {missing}"

    if len(valid_taxa) < 2:
        return False, f"Need at least 2 taxa for MRCA, got {len(valid_taxa)}"

    # Find MRCA
    try:
        mrca = tree.get_common_ancestor(valid_taxa)
    except Exception as e:
        return False, f"Failed to find MRCA: {e}"

    # Check monophyly
    mrca_leaves = set(mrca.get_leaf_names())
    if not set(valid_taxa).issubset(mrca_leaves):
        return False, f"Taxa do not form a monophyletic group"

    # Check if node already has a calibration (starts with 'B(', 'L(', 'U(', or quoted versions)
    if mrca.name and (mrca.name.startswith("'B(") or mrca.name.startswith("'L(") or
                       mrca.name.startswith("'U(") or mrca.name.startswith("B(") or
                       mrca.name.startswith("L(") or mrca.name.startswith("U(")):
        return False, f"Node already has calibration '{mrca.name}'. Skipping [{label}] to avoid invalid Newick."

    # Annotate node with calibration
    # MCMCtree expects calibration in node name
    if mrca.name and mrca.name != "NoName":
        mrca.name = f"{calibration_str}{mrca.name}"
    else:
        mrca.name = calibration_str

    return True, f"Annotated MRCA of {valid_taxa}"


def write_mcmctree_newick(tree: Tree, outfile: Path):
    """
    Write tree in Newick format compatible with MCMCtree.
    MCMCtree expects calibrations as node labels.

    Note: ete3's tree.write() escapes special characters in node names,
    which corrupts calibration format. We need to restore them after writing.
    Also, ete3 doesn't output root node name, so we manually append it.
    """
    # Use format 8 to include internal node names
    newick = tree.write(format=8)

    # ete3 escapes special Newick characters in node names by replacing:
    # '(' -> '_', ')' -> '_', ',' -> '_', ' ' -> '_', ':' -> '_'
    # We need to restore calibration format: 'B_min_max_' -> 'B(min,max)'
    import re

    # Pattern to match calibration in escaped form: 'B_num_num_' or 'L_num_' or 'U_num_'
    # These appear as node names in the tree
    def restore_calibration(match):
        text = match.group(0)
        # Remove surrounding quotes
        inner = text.strip("'")

        if inner.startswith("B_"):
            # Format: B_min_max_ -> B(min,max)
            parts = inner[2:-1].split("_")  # Remove "B_" prefix and trailing "_"
            if len(parts) == 2:
                return f"'B({parts[0]},{parts[1]})'"
        elif inner.startswith("L_"):
            # Format: L_min_ -> L(min)
            parts = inner[2:-1].split("_")
            if len(parts) == 1:
                return f"'L({parts[0]})'"
        elif inner.startswith("U_"):
            # Format: U_max_ -> U(max)
            parts = inner[2:-1].split("_")
            if len(parts) == 1:
                return f"'U({parts[0]})'"

        # If can't parse, return original
        return text

    # Find and replace escaped calibrations
    newick = re.sub(r"'[BLU]_[\d.]+(?:_[\d.]+)?_'", restore_calibration, newick)

    # Handle root node calibration (ete3 doesn't output root name)
    root_name = tree.name
    if root_name:
        # Restore format for root node name too
        root_name_restored = re.sub(r"'[BLU]_[\d.]+(?:_[\d.]+)?_'", restore_calibration, f"'{root_name}'")
        root_name_restored = root_name_restored.strip("'")
        # Remove trailing semicolon, append root name, add semicolon back
        newick = newick.rstrip(';').rstrip() + f"'{root_name_restored}';"

    with open(outfile, 'w') as f:
        f.write(newick)


def main():
    tree_file = Path(snakemake.input.tree)
    calib_file = Path(snakemake.input.calibrations)
    out_tree = Path(snakemake.output.calibrated)
    out_log = Path(snakemake.output.log_file)
    unit = snakemake.params.unit
    
    logs = []
    logs.append(f"Calibration unit: {unit}")
    
    # Load tree
    tree = Tree(str(tree_file))
    logs.append(f"Loaded tree with {len(tree.get_leaves())} taxa")
    
    # Load calibration table
    calib_df = pd.read_csv(calib_file, sep='\t')
    logs.append(f"Loaded {len(calib_df)} calibration constraints")
    
    # Required columns
    required = ['taxa_csv', 'min_age', 'max_age']
    for col in required:
        if col not in calib_df.columns:
            raise ValueError(f"Missing required column: {col}")
    
    success_count = 0
    fail_count = 0
    
    for idx, row in calib_df.iterrows():
        taxa_str = row['taxa_csv']
        min_age = float(row['min_age'])
        max_age = float(row['max_age'])
        label = row.get('label', f"calib_{idx}")
        prior = row.get('prior', None)
        
        # Parse taxa list
        taxa_list = [t.strip() for t in taxa_str.split(',')]
        
        # Format calibration
        calib_str = format_calibration(min_age, max_age, prior)
        
        # Annotate tree
        success, msg = find_mrca_and_annotate(tree, taxa_list, calib_str, label)
        
        if success:
            logs.append(f"OK [{label}]: {msg}, constraint={calib_str}")
            success_count += 1
        else:
            logs.append(f"FAIL [{label}]: {msg}")
            fail_count += 1
    
    logs.append(f"Summary: {success_count} successful, {fail_count} failed")
    
    # Write output
    out_tree.parent.mkdir(parents=True, exist_ok=True)
    write_mcmctree_newick(tree, out_tree)
    
    with open(out_log, 'w') as f:
        f.write('\n'.join(logs))
    
    print(f"Calibrated tree written to {out_tree}")
    print(f"Log written to {out_log}")
    
    if fail_count > 0:
        print(f"WARNING: {fail_count} calibrations failed. Check log for details.")


if __name__ == "__main__":
    main()

