#!/bin/bash
# Setup timetree analysis environment
# Based on README.md - Method 2 (manual environment)
#
# Usage:
#   ./setup_env.sh         # Setup environment
#   ./setup_env.sh clean   # Remove paml directory

set -e

SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

# Handle clean command
if [ "$1" = "clean" ]; then
    echo "Cleaning up..."
    if [ -d "paml" ]; then
        rm -rf paml
        echo "Removed paml directory"
    else
        echo "paml directory not found"
    fi
    exit 0
fi

ENV_NAME="timetree"

echo "=== Setting up TimeTree analysis environment ==="
echo ""

# Check if environment already exists
if conda env list | grep -q "^${ENV_NAME} "; then
    echo "Environment '$ENV_NAME' already exists."
    read -p "Remove and recreate? [y/N] " -n 1 -r
    echo
    if [[ $REPLY =~ ^[Yy]$ ]]; then
        conda env remove -n $ENV_NAME -y
    else
        echo "Skipping environment creation."
        echo "To activate: conda activate $ENV_NAME"
        exit 0
    fi
fi

echo "Creating conda environment: $ENV_NAME"
conda create -n $ENV_NAME python=3.12 -y

echo ""
echo "Installing packages..."

# Install packages
conda run -n $ENV_NAME pip install snakemake biopython ete3 pandas dendropy

conda run -n $ENV_NAME conda install -c conda-forge -c bioconda \
    mafft trimal clipkit seqkit iqtree -y

echo ""
echo "=== Installing IQ2MC modified MCMCtree ==="

# Get environment path
ENV_PATH=$(conda run -n $ENV_NAME printenv CONDA_PREFIX)

# Clone and build modified PAML
if [ -d "paml" ]; then
    echo "paml directory exists, updating..."
    cd paml && git pull && cd ..
else
    echo "Cloning iqtree/paml..."
    git clone https://github.com/iqtree/paml
fi

cd paml/src
echo "Building PAML..."
make clean 2>/dev/null || true
make -j

echo "Installing binaries to $ENV_PATH/bin/"
cp baseml basemlg chi2 codeml evolver infinitesites pamp yn00 mcmctree "$ENV_PATH/bin/"

cd ../..

echo ""
echo "=== Verifying installation ==="
echo ""

conda run -n $ENV_NAME bash -c '
echo "mafft: $(mafft --version 2>&1 | head -1)"
echo "trimal: $(trimal --version 2>&1 | head -1)"
echo "iqtree: $(iqtree --version 2>&1 | head -1)"
echo "mcmctree: $(which mcmctree)"
python -c "import Bio; import ete3; import pandas; print(\"Python packages: OK\")"
'

echo ""
echo "=== Setup complete ==="
echo ""
echo "To activate the environment:"
echo "  conda activate $ENV_NAME"
