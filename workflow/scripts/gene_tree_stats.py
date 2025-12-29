#!/usr/bin/env python3
"""
Generate summary statistics for gene trees.

This script analyzes the gene trees used for coalescent-based species tree inference
and produces summary statistics useful for quality assessment.
"""

from pathlib import Path
import sys


def count_taxa(newick_str: str) -> int:
    """Count number of taxa in a Newick tree string."""
    # Simple approach: count leaf names (between parentheses or commas and colons)
    import re
    # Remove branch lengths
    tree = re.sub(r':[0-9.eE+-]+', '', newick_str)
    # Remove support values in brackets
    tree = re.sub(r'\[[^\]]*\]', '', tree)
    # Count leaves (anything that's not punctuation)
    leaves = re.findall(r'[^(),;:\s]+', tree)
    return len([l for l in leaves if l])


def main():
    gene_trees_file = snakemake.input.gene_trees
    output_file = snakemake.output.stats
    log_file = snakemake.log[0]

    with open(log_file, 'w') as logf:
        logf.write(f"Analyzing gene trees from: {gene_trees_file}\n")

    # Read all gene trees
    with open(gene_trees_file) as f:
        trees = [line.strip() for line in f if line.strip()]

    n_trees = len(trees)
    taxa_counts = [count_taxa(t) for t in trees]

    # Calculate statistics
    min_taxa = min(taxa_counts) if taxa_counts else 0
    max_taxa = max(taxa_counts) if taxa_counts else 0
    avg_taxa = sum(taxa_counts) / len(taxa_counts) if taxa_counts else 0

    # Write statistics
    with open(output_file, 'w') as f:
        f.write("Gene Tree Summary Statistics\n")
        f.write("=" * 40 + "\n\n")
        f.write(f"Total number of gene trees: {n_trees}\n")
        f.write(f"Taxa per tree:\n")
        f.write(f"  - Minimum: {min_taxa}\n")
        f.write(f"  - Maximum: {max_taxa}\n")
        f.write(f"  - Average: {avg_taxa:.1f}\n")
        f.write("\n")
        f.write("Taxa count distribution:\n")

        # Distribution histogram
        from collections import Counter
        dist = Counter(taxa_counts)
        for taxa, count in sorted(dist.items()):
            bar = "#" * min(count, 50)
            f.write(f"  {taxa:3d} taxa: {count:4d} trees {bar}\n")

    with open(log_file, 'a') as logf:
        logf.write(f"Wrote statistics to: {output_file}\n")


if __name__ == "__main__":
    main()
