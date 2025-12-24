#!/usr/bin/env python3
"""
Post-process MCMCtree output to generate final time tree files.
"""

from pathlib import Path
import re


def parse_figtree(figtree_path: Path) -> str:
    """
    Parse FigTree.tre from MCMCtree and extract tree with node ages.
    """
    with open(figtree_path) as f:
        content = f.read()
    
    # FigTree format contains tree block
    # Extract the tree string
    tree_match = re.search(r'tree\s+\w+\s*=\s*\[&[^\]]*\]\s*(.+);', content, re.IGNORECASE)
    if tree_match:
        return tree_match.group(1) + ';'
    
    # Try simpler format
    for line in content.strip().split('\n'):
        line = line.strip()
        if line.startswith('(') and line.endswith(';'):
            return line
    
    return content.strip()


def convert_to_nexus(newick: str, title: str = "TimeTree") -> str:
    """Convert Newick to NEXUS format for FigTree."""
    return f"""#NEXUS
BEGIN TREES;
    TREE {title} = {newick}
END;
"""


def main():
    figtree_in = Path(snakemake.input.figtree)
    mcmc_in = Path(snakemake.input.mcmc)
    out_nwk = Path(snakemake.output.nwk)
    out_nex = Path(snakemake.output.nex)
    
    out_nwk.parent.mkdir(parents=True, exist_ok=True)
    
    # Parse FigTree output
    tree_str = parse_figtree(figtree_in)
    
    # Write Newick
    with open(out_nwk, 'w') as f:
        f.write(tree_str + '\n')
    
    # Write NEXUS
    nexus_content = convert_to_nexus(tree_str)
    with open(out_nex, 'w') as f:
        f.write(nexus_content)
    
    # Report MCMC summary if available
    if mcmc_in.exists():
        with open(mcmc_in) as f:
            lines = f.readlines()
        n_samples = len([l for l in lines if not l.startswith('#') and l.strip()])
        print(f"MCMC samples: {n_samples}")
    
    print(f"Final time tree written to:")
    print(f"  Newick: {out_nwk}")
    print(f"  NEXUS:  {out_nex}")


if __name__ == "__main__":
    main()

