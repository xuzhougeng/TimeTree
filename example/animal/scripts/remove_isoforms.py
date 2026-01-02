#!/usr/bin/env python3
"""
Remove isoforms from Ensembl FASTA files.
For each gene, keep only the longest transcript.
"""

import gzip
import sys
import re
from pathlib import Path
from collections import defaultdict


def parse_gene_id(header):
    """Extract gene ID from Ensembl FASTA header."""
    match = re.search(r'gene:(\S+)', header)
    if match:
        return match.group(1)
    # Fallback: try to extract gene ID from sequence ID (e.g., Vitvi12g02389_P001 -> Vitvi12g02389)
    seq_id = header.split()[0].lstrip('>')
    match = re.match(r'^([^_]+)', seq_id)
    if match:
        return match.group(1)
    return None


def read_fasta(file_path):
    """Read FASTA file (gzipped or plain text)."""
    opener = gzip.open if file_path.endswith('.gz') else open
    sequences = []

    with opener(file_path, 'rt') as f:
        current_header = None
        current_seq = []

        for line in f:
            line = line.rstrip()
            if line.startswith('>'):
                if current_header:
                    sequences.append((current_header, ''.join(current_seq)))
                current_header = line
                current_seq = []
            else:
                current_seq.append(line)

        if current_header:
            sequences.append((current_header, ''.join(current_seq)))

    return sequences


def remove_isoforms(sequences):
    """Keep only the longest transcript for each gene."""
    gene_to_transcripts = defaultdict(list)

    # Group sequences by gene
    for header, seq in sequences:
        gene_id = parse_gene_id(header)
        if gene_id:
            gene_to_transcripts[gene_id].append((header, seq, len(seq)))
        else:
            print(f"Warning: Could not parse gene ID from header: {header[:80]}", file=sys.stderr)
            gene_to_transcripts[header].append((header, seq, len(seq)))

    # Select longest transcript for each gene
    longest_transcripts = []
    for gene_id, transcripts in gene_to_transcripts.items():
        longest = max(transcripts, key=lambda x: x[2])
        longest_transcripts.append((longest[0], longest[1]))

    return longest_transcripts


def write_fasta(sequences, output_path):
    """Write sequences to FASTA file (plain text, optimized for speed)."""
    # Always write as plain text for better performance
    with open(output_path, 'w') as f:
        for header, seq in sequences:
            f.write(f"{header}\n")
            # Write sequence in 60-character lines
            for i in range(0, len(seq), 60):
                f.write(f"{seq[i:i+60]}\n")


def main():
    if len(sys.argv) != 3:
        print("Usage: remove_isoforms.py <input.fa[.gz]> <output.fa[.gz]>", file=sys.stderr)
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]

    print(f"Reading sequences from {input_file}...")
    sequences = read_fasta(input_file)
    print(f"  Found {len(sequences)} sequences")

    print("Removing isoforms (keeping longest transcript per gene)...")
    longest = remove_isoforms(sequences)
    print(f"  Kept {len(longest)} sequences")
    print(f"  Removed {len(sequences) - len(longest)} isoforms")

    print(f"Writing to {output_file}...")
    write_fasta(longest, output_file)
    print("Done!")


if __name__ == '__main__':
    main()
