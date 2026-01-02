#!/bin/bash
# Remove isoforms from already downloaded CDS and PEP sequences
# This script is useful if you already have downloaded data and want to remove isoforms

set -e

SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

REMOVE_ISOFORMS_SCRIPT="$SCRIPT_DIR/remove_isoforms.py"

if [ ! -f "$REMOVE_ISOFORMS_SCRIPT" ]; then
    echo "Error: remove_isoforms.py not found at $REMOVE_ISOFORMS_SCRIPT"
    exit 1
fi

# Check if data directories exist
if [ ! -d "data/cds" ] || [ ! -d "data/pep" ]; then
    echo "Error: data/cds or data/pep directory not found"
    echo "Please run ./scripts/download_plant_data.sh first"
    exit 1
fi

echo "=== Removing isoforms from existing data ==="
echo ""

# Create temporary directory for processed files
mkdir -p data/cds_filtered data/pep_filtered

# Process all CDS files
echo "Processing CDS files..."
for cds_file in data/cds/*.fa.gz; do
    if [ -f "$cds_file" ]; then
        basename=$(basename "$cds_file" .fa.gz)
        echo "  Processing $basename..."
        python3 "$REMOVE_ISOFORMS_SCRIPT" "$cds_file" "data/cds_filtered/${basename}.fa"
    fi
done
echo ""

# Process all PEP files
echo "Processing PEP files..."
for pep_file in data/pep/*.fa.gz; do
    if [ -f "$pep_file" ]; then
        basename=$(basename "$pep_file" .fa.gz)
        echo "  Processing $basename..."
        python3 "$REMOVE_ISOFORMS_SCRIPT" "$pep_file" "data/pep_filtered/${basename}.fa"
    fi
done
echo ""

# Backup original files and replace with filtered ones
echo "Backing up original files and replacing with filtered versions..."
if [ -d "data/cds_original" ]; then
    echo "Warning: data/cds_original already exists, removing it..."
    rm -rf data/cds_original
fi
if [ -d "data/pep_original" ]; then
    echo "Warning: data/pep_original already exists, removing it..."
    rm -rf data/pep_original
fi

mv data/cds data/cds_original
mv data/cds_filtered data/cds
mv data/pep data/pep_original
mv data/pep_filtered data/pep

echo ""
echo "=== Isoform removal complete ==="
echo "Original files backed up to:"
echo "  - data/cds_original/"
echo "  - data/pep_original/"
echo ""
echo "Filtered files are now in:"
echo "  - data/cds/"
echo "  - data/pep/"
echo ""
echo "Summary:"
echo "CDS files:"
ls -lh data/cds/
echo ""
echo "PEP files:"
ls -lh data/pep/
