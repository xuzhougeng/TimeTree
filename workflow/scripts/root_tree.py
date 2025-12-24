#!/usr/bin/env python3
"""
Root a phylogenetic tree using outgroup or midpoint rooting.
"""

from pathlib import Path
from ete3 import Tree


def main():
    tree_file = Path(snakemake.input.tree)
    out_file = Path(snakemake.output.rooted)
    outgroup = snakemake.params.outgroup
    
    # Load tree
    tree = Tree(str(tree_file))
    
    # Get leaf names for validation
    leaf_names = set(tree.get_leaf_names())
    
    if outgroup and outgroup != []:
        # Filter outgroup to only include taxa present in tree
        valid_outgroup = [t for t in outgroup if t in leaf_names]
        
        if valid_outgroup:
            if len(valid_outgroup) == 1:
                tree.set_outgroup(valid_outgroup[0])
            else:
                # Get MRCA of outgroup taxa and root there
                ancestor = tree.get_common_ancestor(valid_outgroup)
                tree.set_outgroup(ancestor)
            print(f"Rooted tree using outgroup: {valid_outgroup}")
        else:
            print(f"WARNING: No valid outgroup taxa found. Using midpoint rooting.")
            midpoint = tree.get_midpoint_outgroup()
            if midpoint:
                tree.set_outgroup(midpoint)
    else:
        # Midpoint rooting
        midpoint = tree.get_midpoint_outgroup()
        if midpoint:
            tree.set_outgroup(midpoint)
        print("Rooted tree using midpoint rooting")
    
    # Write rooted tree
    out_file.parent.mkdir(parents=True, exist_ok=True)
    tree.write(outfile=str(out_file), format=1)
    
    print(f"Rooted tree written to {out_file}")


if __name__ == "__main__":
    main()

