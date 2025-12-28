#!/bin/bash
# Download CDS and PEP sequences for 5 plant species from Ensembl Plants
# Species: Arabidopsis thaliana, Oryza sativa, Zea mays, Solanum lycopersicum, Vitis vinifera

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
)

echo "=== Downloading CDS and PEP sequences for 5 plant species ==="
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
echo "To decompress files, run:"
echo "  gunzip data/cds/*.gz data/pep/*.gz"
