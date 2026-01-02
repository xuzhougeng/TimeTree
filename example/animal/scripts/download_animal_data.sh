#!/bin/bash
# Download CDS and PEP sequences for 8 animal species from Ensembl
# Species: Human, Mouse, Cow, Dog, Chicken, Zebrafish, Fruit fly, C. elegans
# Vertebrates from Ensembl main, Invertebrates from Ensembl Metazoa

set -e

# Create output directories
SCRIPT_DIR=$(dirname "$0")
WORKDIR=$(dirname "$SCRIPT_DIR")
cd "$WORKDIR"

mkdir -p data/cds data/pep

# Ensembl FTP base URLs
ENSEMBL_URL="https://ftp.ensembl.org/pub/release-113/fasta"
METAZOA_URL="https://ftp.ensemblgenomes.ebi.ac.uk/pub/metazoa/release-60/fasta"

# Define vertebrate species: species_name, dir_name, file_prefix
VERTEBRATES=(
    "Homo_sapiens:homo_sapiens:Homo_sapiens.GRCh38"
    "Mus_musculus:mus_musculus:Mus_musculus.GRCm39"
    "Bos_taurus:bos_taurus:Bos_taurus.ARS-UCD1.3"
    "Canis_lupus_familiaris:canis_lupus_familiaris:Canis_lupus_familiaris.ROS_Cfam_1.0"
    "Gallus_gallus:gallus_gallus:Gallus_gallus.bGalGal1.mat.broiler.GRCg7b"
    "Danio_rerio:danio_rerio:Danio_rerio.GRCz11"
)

# Define invertebrate species (from Ensembl Metazoa)
INVERTEBRATES=(
    "Drosophila_melanogaster:drosophila_melanogaster:Drosophila_melanogaster.BDGP6.46"
    "Caenorhabditis_elegans:caenorhabditis_elegans:Caenorhabditis_elegans.WBcel235"
)

echo "=== Downloading CDS and PEP sequences for 8 animal species ==="
echo "Data source: Ensembl Release 113 (Vertebrates), Ensembl Metazoa Release 60 (Invertebrates)"
echo ""

# Download vertebrates from Ensembl main
echo "--- Downloading Vertebrates from Ensembl ---"
for entry in "${VERTEBRATES[@]}"; do
    IFS=':' read -r species dir_name file_prefix <<< "$entry"

    echo "Processing: $species"

    # Download CDS
    cds_out="data/cds/${species}.cds.fa.gz"
    if [ -f "$cds_out" ] && [ -s "$cds_out" ]; then
        echo "  CDS: already exists, skipping"
    else
        cds_file="${file_prefix}.cds.all.fa.gz"
        cds_url="${ENSEMBL_URL}/${dir_name}/cds/${cds_file}"
        echo "  Downloading CDS: $cds_file"
        wget -q --show-progress -O "$cds_out" "$cds_url"
    fi

    # Download PEP
    pep_out="data/pep/${species}.pep.fa.gz"
    if [ -f "$pep_out" ] && [ -s "$pep_out" ]; then
        echo "  PEP: already exists, skipping"
    else
        pep_file="${file_prefix}.pep.all.fa.gz"
        pep_url="${ENSEMBL_URL}/${dir_name}/pep/${pep_file}"
        echo "  Downloading PEP: $pep_file"
        wget -q --show-progress -O "$pep_out" "$pep_url"
    fi

    echo "  Done!"
    echo ""
done

# Download invertebrates from Ensembl Metazoa
echo "--- Downloading Invertebrates from Ensembl Metazoa ---"
for entry in "${INVERTEBRATES[@]}"; do
    IFS=':' read -r species dir_name file_prefix <<< "$entry"

    echo "Processing: $species"

    # Download CDS
    cds_out="data/cds/${species}.cds.fa.gz"
    if [ -f "$cds_out" ] && [ -s "$cds_out" ]; then
        echo "  CDS: already exists, skipping"
    else
        cds_file="${file_prefix}.cds.all.fa.gz"
        cds_url="${METAZOA_URL}/${dir_name}/cds/${cds_file}"
        echo "  Downloading CDS: $cds_file"
        wget -q --show-progress -O "$cds_out" "$cds_url"
    fi

    # Download PEP
    pep_out="data/pep/${species}.pep.fa.gz"
    if [ -f "$pep_out" ] && [ -s "$pep_out" ]; then
        echo "  PEP: already exists, skipping"
    else
        pep_file="${file_prefix}.pep.all.fa.gz"
        pep_url="${METAZOA_URL}/${dir_name}/pep/${pep_file}"
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

ALL_SPECIES=("${VERTEBRATES[@]}" "${INVERTEBRATES[@]}")

for entry in "${ALL_SPECIES[@]}"; do
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
