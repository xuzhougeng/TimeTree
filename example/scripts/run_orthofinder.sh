#!/bin/bash
# Run OrthoFinder on the 5 plant species protein sequences
# Requires: orthofinder (conda install -c bioconda orthofinder)

set -e

SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

# Check if pep files exist
if [ ! -d "data/pep" ] || [ -z "$(ls -A data/pep/*.gz 2>/dev/null)" ]; then
    echo "Error: No protein files found in data/pep/"
    echo "Please run scripts/download_plant_data.sh first"
    exit 1
fi

# Create orthofinder input directory
mkdir -p orthofinder_input

echo "=== Preparing protein sequences for OrthoFinder ==="

# Decompress and prepare protein files
for f in data/pep/*.pep.fa.gz; do
    species=$(basename "$f" .pep.fa.gz)
    echo "Processing: $species"

    # Decompress and clean headers (keep only gene ID)
    zcat "$f" | sed 's/\.[0-9]*$//' > "orthofinder_input/${species}.fa"
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
