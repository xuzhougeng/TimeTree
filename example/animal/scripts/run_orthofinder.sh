#!/bin/bash
# Run OrthoFinder on the 8 animal species protein sequences
# Requires: orthofinder (conda install -c bioconda orthofinder)

set -e

SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

# Check if pep files exist
if [ ! -d "data/pep" ]; then
    echo "Error: data/pep/ directory not found"
    echo "Please run scripts/download_animal_data.sh first"
    exit 1
fi

# Check for both .fa and .fa.gz files
pep_files=$(ls data/pep/*.pep.fa 2>/dev/null || ls data/pep/*.pep.fa.gz 2>/dev/null || true)
if [ -z "$pep_files" ]; then
    echo "Error: No protein files found in data/pep/"
    echo "Please run scripts/download_animal_data.sh first"
    exit 1
fi

# Create orthofinder input directory
mkdir -p orthofinder_input

echo "=== Preparing protein sequences for OrthoFinder ==="

# Process protein files (both .fa and .fa.gz)
for f in data/pep/*.pep.fa data/pep/*.pep.fa.gz; do
    # Skip if file doesn't exist (in case only one format is present)
    [ -e "$f" ] || continue

    # Extract species name
    if [[ "$f" == *.fa.gz ]]; then
        species=$(basename "$f" .pep.fa.gz)
    else
        species=$(basename "$f" .pep.fa)
    fi

    echo "Processing: $species"

    # Handle both compressed and uncompressed files
    if [[ "$f" == *.gz ]]; then
        # Decompress and clean headers
        zcat "$f" | sed 's/\.[0-9]*$//' > "orthofinder_input/${species}.fa"
    else
        # Just clean headers (no decompression needed)
        cat "$f" | sed 's/\.[0-9]*$//' > "orthofinder_input/${species}.fa"
    fi
done

# Check if output directory exists
if [ -d "orthofinder_analysis" ]; then
    echo ""
    echo "Warning: orthofinder_analysis/ already exists"
    read -p "Remove and rerun? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        rm -rf orthofinder_analysis
    else
        echo "Aborted."
        exit 0
    fi
fi

echo ""
echo "=== Running OrthoFinder ==="
echo "Input directory: orthofinder_input/"
echo ""

# Run OrthoFinder
# -f: input directory
# -t: number of threads
# -a: number of threads for analysis
# -S: sequence search program (diamond is faster)
conda run -n of orthofinder -f orthofinder_input \
    -t 64 \
    -a 4 \
    -S diamond \
    -o orthofinder_analysis

echo ""
echo "=== OrthoFinder Complete ==="
echo ""
echo "Results directory:"
ls -d orthofinder_analysis/Results_*

echo ""
echo "Single-copy orthologs:"
cat orthofinder_analysis/Results_*/Orthogroups/Orthogroups_SingleCopyOrthologues.txt | wc -l
