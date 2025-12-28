#!/bin/bash
# Download CDS and PEP sequences for 8 plant species from Ensembl Plants
# Species: Arabidopsis thaliana, Oryza sativa, Zea mays, Solanum lycopersicum,
#          Vitis vinifera, Marchantia polymorpha, Brassica rapa, Amborella trichopoda

set -e

# Create output directories
SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

mkdir -p data/cds data/pep

# Ensembl Plants FTP base URL (release 60)
BASE_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta"

# Define species info: species_name, dir_name, file_prefix
SPECIES_LIST=(
    "Arabidopsis_thaliana:arabidopsis_thaliana:Arabidopsis_thaliana.TAIR10"
    "Oryza_sativa:oryza_sativa:Oryza_sativa.IRGSP-1.0"
    "Zea_mays:zea_mays:Zea_mays.Zm-B73-REFERENCE-NAM-5.0"
    "Solanum_lycopersicum:solanum_lycopersicum:Solanum_lycopersicum.SL3.0"
    "Vitis_vinifera:vitis_vinifera:Vitis_vinifera.PN40024.v4"
    "Marchantia_polymorpha:marchantia_polymorpha:Marchantia_polymorpha.Marchanta_polymorpha_v1"
    "Brassica_rapa:brassica_rapa:Brassica_rapa.Brapa_1.0"
    "Amborella_trichopoda:amborella_trichopoda:Amborella_trichopoda.AMTR1.0"
)

echo "=== Downloading CDS and PEP sequences for 8 plant species ==="
echo "Data source: Ensembl Plants Release 60"
echo ""

for entry in "${SPECIES_LIST[@]}"; do
    IFS=':' read -r species dir_name file_prefix <<< "$entry"

    echo "Processing: $species"

    # Download CDS
    cds_out="data/cds/${species}.cds.fa.gz"
    if [ -f "$cds_out" ] && [ -s "$cds_out" ]; then
        echo "  CDS: already exists, skipping"
    else
        cds_file="${file_prefix}.cds.all.fa.gz"
        cds_url="${BASE_URL}/${dir_name}/cds/${cds_file}"
        echo "  Downloading CDS: $cds_file"
        wget -q --show-progress -O "$cds_out" "$cds_url"
    fi

    # Download PEP
    pep_out="data/pep/${species}.pep.fa.gz"
    if [ -f "$pep_out" ] && [ -s "$pep_out" ]; then
        echo "  PEP: already exists, skipping"
    else
        pep_file="${file_prefix}.pep.all.fa.gz"
        pep_url="${BASE_URL}/${dir_name}/pep/${pep_file}"
        echo "  Downloading PEP: $pep_file"
        wget -q --show-progress -O "$pep_out" "$pep_url"
    fi

    echo "  Done!"
    echo ""
done

echo "=== Download complete ==="
echo ""
echo "CDS files:"
ls -lh data/cds/
echo ""
echo "PEP files:"
ls -lh data/pep/
echo ""

# Remove isoforms (keep only longest transcript per gene)
echo "=== Removing isoforms (keeping longest transcript per gene) ==="
echo ""

REMOVE_ISOFORMS_SCRIPT="$SCRIPT_DIR/remove_isoforms.py"

if [ ! -f "$REMOVE_ISOFORMS_SCRIPT" ]; then
    echo "Error: remove_isoforms.py not found at $REMOVE_ISOFORMS_SCRIPT"
    exit 1
fi

# Create temporary directory for processed files
mkdir -p data/cds_filtered data/pep_filtered

for entry in "${SPECIES_LIST[@]}"; do
    IFS=':' read -r species dir_name file_prefix <<< "$entry"

    echo "Processing $species..."

    # Remove isoforms from CDS (output uncompressed for speed)
    cds_in="data/cds/${species}.cds.fa.gz"
    cds_out="data/cds_filtered/${species}.cds.fa"
    if [ -f "$cds_in" ]; then
        echo "  Removing CDS isoforms..."
        python3 "$REMOVE_ISOFORMS_SCRIPT" "$cds_in" "$cds_out"
    fi

    # Remove isoforms from PEP (output uncompressed for speed)
    pep_in="data/pep/${species}.pep.fa.gz"
    pep_out="data/pep_filtered/${species}.pep.fa"
    if [ -f "$pep_in" ]; then
        echo "  Removing PEP isoforms..."
        python3 "$REMOVE_ISOFORMS_SCRIPT" "$pep_in" "$pep_out"
    fi

    echo ""
done

# Replace original files with filtered ones
echo "Replacing original files with filtered versions..."
mv data/cds data/cds_original
mv data/cds_filtered data/cds
mv data/pep data/pep_original
mv data/pep_filtered data/pep

echo ""
echo "=== Isoform removal complete ==="
echo "Original files (compressed) backed up to data/cds_original/ and data/pep_original/"
echo "Filtered files (uncompressed) are now in data/cds/ and data/pep/"
echo ""
echo "Filtered CDS files:"
ls -lh data/cds/
echo ""
echo "Filtered PEP files:"
ls -lh data/pep/
