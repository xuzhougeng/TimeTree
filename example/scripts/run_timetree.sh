#!/bin/bash
# Run TimeTree workflow
# Requires: timetree conda environment (run setup_env.sh first)

set -e

# Get absolute paths
SCRIPT_DIR=$(cd "$(dirname "$0")" && pwd)
WORKDIR=$(dirname "$SCRIPT_DIR")
TIMETREE_DIR=$(dirname "$WORKDIR")

cd "$WORKDIR"

echo "=== TimeTree Workflow ==="
echo "Working directory: $WORKDIR"
echo "Workflow directory: $TIMETREE_DIR"
echo ""

# Check prerequisites
echo "Checking prerequisites..."

# Check OrthoFinder results
RESULTS_DIR=$(ls -d orthofinder_analysis/Results_* 2>/dev/null | head -1)
if [ -z "$RESULTS_DIR" ]; then
    echo "Error: No OrthoFinder results found in orthofinder_analysis/"
    echo "Please run scripts/run_orthofinder.sh first"
    exit 1
fi
echo "  OrthoFinder results: $RESULTS_DIR"

# Check calibrations
if [ ! -f "calibrations.tsv" ]; then
    echo "Error: calibrations.tsv not found"
    exit 1
fi
echo "  Calibrations: calibrations.tsv"

# Check config
if [ ! -f "config.yaml" ]; then
    echo "Error: config.yaml not found"
    exit 1
fi
echo "  Config: config.yaml"

# Update results_dir in config.yaml if needed
CURRENT_RESULTS=$(grep "results_dir:" config.yaml | sed 's/.*: *"\(.*\)"/\1/')
if [ "$CURRENT_RESULTS" != "$RESULTS_DIR" ]; then
    echo ""
    echo "Updating results_dir in config.yaml to: $RESULTS_DIR"
    sed -i "s|results_dir:.*|results_dir: \"$RESULTS_DIR\"|" config.yaml
fi

echo ""
echo "=== Running Snakemake workflow ==="
echo ""

# Run snakemake from the TimeTree directory with example config
conda run -n timetree snakemake \
    --snakefile "$TIMETREE_DIR/Snakefile" \
    --configfile config.yaml \
    --cores 64 \
    --rerun-incomplete \
    "$@"

echo ""
echo "=== Workflow complete ==="
echo ""
echo "Results in: results/timetree/"
ls -la results/timetree/ 2>/dev/null || echo "(results directory not yet created)"
