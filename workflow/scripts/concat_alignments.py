#!/usr/bin/env python3
"""
Concatenate multiple sequence alignments into a supermatrix.
Generate partition file in NEXUS format for IQ-TREE.
"""

from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


def parse_alignment(aln_file: Path) -> dict[str, str]:
    """Parse FASTA alignment into dict of species -> sequence."""
    seqs = {}
    for record in SeqIO.parse(aln_file, "fasta"):
        seqs[record.id] = str(record.seq)
    return seqs


def main():
    aln_files = [Path(f) for f in snakemake.input.alignments]
    out_phy = Path(snakemake.output.supermatrix)
    out_part = Path(snakemake.output.partitions)
    out_fasta = Path(snakemake.output.fasta)
    
    out_phy.parent.mkdir(parents=True, exist_ok=True)
    
    # Collect all species and alignments
    all_species = set()
    alignments = {}
    
    for aln_file in sorted(aln_files):
        og_id = aln_file.stem.replace('.aln', '')
        seqs = parse_alignment(aln_file)
        alignments[og_id] = seqs
        all_species.update(seqs.keys())
    
    all_species = sorted(all_species)
    
    # Build supermatrix
    supermatrix = defaultdict(str)
    partitions = []
    pos = 1
    
    for og_id in sorted(alignments.keys()):
        seqs = alignments[og_id]
        # Get alignment length from first sequence
        aln_len = len(next(iter(seqs.values())))
        
        for sp in all_species:
            if sp in seqs:
                supermatrix[sp] += seqs[sp]
            else:
                # Missing species: fill with gaps
                supermatrix[sp] += '-' * aln_len
        
        # Record partition
        partitions.append((og_id, pos, pos + aln_len - 1))
        pos += aln_len
    
    total_len = len(next(iter(supermatrix.values())))
    n_taxa = len(all_species)
    
    # Write PHYLIP format
    with open(out_phy, 'w') as f:
        f.write(f"{n_taxa} {total_len}\n")
        for sp in all_species:
            # PHYLIP format: name padded to 10 chars or separated by space
            f.write(f"{sp}  {supermatrix[sp]}\n")
    
    # Write FASTA format
    with open(out_fasta, 'w') as f:
        for sp in all_species:
            f.write(f">{sp}\n{supermatrix[sp]}\n")
    
    # Write partition file in NEXUS format
    with open(out_part, 'w') as f:
        f.write("#nexus\n")
        f.write("begin sets;\n")
        for og_id, start, end in partitions:
            f.write(f"  charset {og_id} = {start}-{end};\n")
        f.write("end;\n")
    
    print(f"Supermatrix: {n_taxa} taxa, {total_len} sites, {len(partitions)} partitions")


if __name__ == "__main__":
    main()

