#!/usr/bin/env python3
"""
Extract single-copy orthogroup (SCO) FASTA files from OrthoFinder results.
If Orthogroup_Sequences exists, use those; otherwise extract from proteomes.
"""

import sys
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO


def parse_sco_list(sco_file: Path) -> list[str]:
    """Parse Orthogroups_SingleCopyOrthologues.txt to get OG IDs."""
    ogs = []
    with open(sco_file) as f:
        for line in f:
            og = line.strip()
            if og:
                ogs.append(og)
    return ogs


def parse_orthogroups_tsv(results_dir: Path) -> dict[str, dict[str, str]]:
    """Parse Orthogroups.tsv to get gene IDs per species per OG."""
    og_file = results_dir / "Orthogroups" / "Orthogroups.tsv"
    og_data = {}
    
    with open(og_file) as f:
        header = f.readline().strip().split('\t')
        species = header[1:]  # First column is Orthogroup
        
        for line in f:
            parts = line.strip().split('\t')
            og_id = parts[0]
            og_data[og_id] = {}
            
            for i, sp in enumerate(species):
                genes = parts[i + 1] if i + 1 < len(parts) else ""
                if genes:
                    # Genes are comma-separated
                    og_data[og_id][sp] = [g.strip() for g in genes.split(',')]
                else:
                    og_data[og_id][sp] = []
    
    return og_data, species


def load_proteomes(proteomes_dir: Path) -> dict[str, dict[str, str]]:
    """Load all proteome sequences indexed by gene ID."""
    proteomes = {}
    for fasta_file in proteomes_dir.glob("*.fa*"):
        species = fasta_file.stem
        proteomes[species] = {}
        for record in SeqIO.parse(fasta_file, "fasta"):
            proteomes[species][record.id] = str(record.seq)
    return proteomes


def main():
    # Snakemake provides input/output/params/log
    sco_list_file = Path(snakemake.input.sco_list)
    og_sequences_dir = Path(snakemake.input.og_dir)
    proteomes_dir = Path(snakemake.input.proteomes)
    output_dir = Path(snakemake.output[0])
    log_file = Path(snakemake.output.log_file)
    min_taxa = snakemake.params.min_taxa
    allow_missing = snakemake.params.allow_missing
    
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Parse SCO list
    sco_ids = parse_sco_list(sco_list_file)
    
    # Get results directory from sco_list path
    results_dir = sco_list_file.parent.parent
    
    logs = []
    logs.append(f"Found {len(sco_ids)} single-copy orthogroups")
    
    extracted = 0
    skipped = 0
    
    # Check if Orthogroup_Sequences exists with FASTAs
    use_og_sequences = og_sequences_dir.exists() and list(og_sequences_dir.glob("OG*.fa"))
    
    if use_og_sequences:
        logs.append(f"Using sequences from {og_sequences_dir}")
        for og_id in sco_ids:
            og_fasta = og_sequences_dir / f"{og_id}.fa"
            if og_fasta.exists():
                # Read and check taxa count
                records = list(SeqIO.parse(og_fasta, "fasta"))
                if len(records) >= min_taxa or allow_missing:
                    out_file = output_dir / f"{og_id}.faa"
                    # Rename sequences to species name (remove gene ID suffix)
                    with open(out_file, 'w') as out:
                        for rec in records:
                            # OrthoFinder names: Species_GeneID or just GeneID
                            # Try to extract species from sequence name
                            seq_name = rec.id.split('_')[0] if '_' in rec.id else rec.id
                            out.write(f">{seq_name}\n{rec.seq}\n")
                    extracted += 1
                else:
                    logs.append(f"SKIP {og_id}: only {len(records)} taxa (min={min_taxa})")
                    skipped += 1
            else:
                logs.append(f"SKIP {og_id}: FASTA not found")
                skipped += 1
    else:
        logs.append(f"Orthogroup_Sequences not found, extracting from proteomes")
        og_data, species_list = parse_orthogroups_tsv(results_dir)
        proteomes = load_proteomes(proteomes_dir)
        logs.append(f"Loaded {len(proteomes)} proteomes with species: {list(proteomes.keys())}")
        
        for og_id in sco_ids:
            if og_id not in og_data:
                logs.append(f"SKIP {og_id}: not in Orthogroups.tsv")
                skipped += 1
                continue
            
            og_genes = og_data[og_id]
            sequences = []
            
            for sp in species_list:
                genes = og_genes.get(sp, [])
                if len(genes) == 1:
                    gene_id = genes[0]
                    if sp in proteomes and gene_id in proteomes[sp]:
                        sequences.append((sp, proteomes[sp][gene_id]))
                    else:
                        logs.append(f"WARN {og_id}: gene {gene_id} not found in {sp} proteome")
            
            if len(sequences) >= min_taxa or allow_missing:
                out_file = output_dir / f"{og_id}.faa"
                with open(out_file, 'w') as out:
                    for sp, seq in sequences:
                        out.write(f">{sp}\n{seq}\n")
                extracted += 1
            else:
                logs.append(f"SKIP {og_id}: only {len(sequences)} taxa (min={min_taxa})")
                skipped += 1
    
    logs.append(f"Extracted: {extracted}, Skipped: {skipped}")
    
    with open(log_file, 'w') as f:
        f.write('\n'.join(logs))
    
    print(f"Extracted {extracted} SCO FASTAs to {output_dir}")


if __name__ == "__main__":
    main()

