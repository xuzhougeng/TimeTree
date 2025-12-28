# TimeTree Example: 8 Plant Species

This example demonstrates the TimeTree workflow using 8 plant species from Ensembl Plants, representing major plant lineages.

## Species

| Species | Common Name | Lineage |
|---------|-------------|---------|
| Marchantia polymorpha | Liverwort | Bryophyte (non-vascular) |
| Amborella trichopoda | Amborella | Basal angiosperm |
| Arabidopsis thaliana | Thale cress | Eudicot (Brassicaceae) |
| Brassica rapa | Chinese cabbage | Eudicot (Brassicaceae) |
| Populus trichocarpa | Poplar | Eudicot (Salicaceae) |
| Solanum lycopersicum | Tomato | Eudicot (Solanaceae) |
| Oryza sativa | Rice | Monocot (Poaceae) |
| Zea mays | Maize | Monocot (Poaceae) |

## Fossil Calibrations

| Calibration Point | Taxa | Age (Mya) | Reference |
|-------------------|------|-----------|-----------|
| Monocot-Dicot Split | Arabidopsis, Oryza | 142-163 | Magallon et al. (2015) |
| Poaceae Crown | Oryza, Zea | 41-52 | Prasad et al. (2011) |

## Directory Structure

```
example/
├── calibrations.tsv          # Fossil calibration data
├── config.yaml               # Workflow configuration
├── data/
│   ├── cds/                  # CDS sequences
│   └── pep/                  # Protein sequences
├── orthofinder_input/        # OrthoFinder input (prepared)
├── orthofinder_analysis/     # OrthoFinder results
├── results/                  # TimeTree output
└── scripts/
    ├── download_plant_data.sh
    ├── remove_isoforms.py        # Remove transcript isoforms
    ├── run_orthofinder.sh
    ├── setup_env.sh
    └── run_timetree.sh
```

## Quick Start

### Step 1: Download Data

```bash
./scripts/download_plant_data.sh
```

Downloads CDS and protein sequences from Ensembl Plants Release 60, then automatically removes isoforms (keeping only the longest transcript per gene).

**What happens:**
1. Downloads `.cds.all.fa.gz` and `.pep.all.fa.gz` files from Ensembl Plants
2. Removes isoforms using `remove_isoforms.py` - keeps only the longest transcript for each gene
3. Outputs uncompressed `.fa` files for faster downstream processing
4. Backs up original compressed files to `data/cds_original/` and `data/pep_original/`
5. Filtered files are saved as uncompressed `.fa` in `data/cds/` and `data/pep/`

**Why remove isoforms?**
Ensembl sequences include all transcript isoforms (alternative splicing variants) of each gene. For phylogenetic analysis, we need only one representative sequence per gene to avoid redundancy and ensure accurate ortholog detection.

**Alternative: Filter existing data**
If you already have downloaded data without isoform removal, you can run:
```bash
./scripts/filter_existing_data.sh
```

### Step 2: Run OrthoFinder

```bash
# conda create -n of -c conda-forge -c bioconda orthofinder
./scripts/run_orthofinder.sh
```

Identifies orthogroups and single-copy orthologs. Requires `of` conda environment with OrthoFinder installed.

### Step 3: Setup Environment

```bash
./scripts/setup_env.sh
```

Creates `timetree` conda environment with:
- Snakemake
- MAFFT, trimAl, IQ-TREE
- Modified MCMCtree (from iqtree/paml)

### Step 4: Run TimeTree

```bash
./scripts/run_timetree.sh
```

Executes the complete workflow:
1. Extract single-copy orthologs
2. Multiple sequence alignment (MAFFT)
3. Alignment trimming (trimAl)
4. Concatenate to supermatrix
5. ML tree inference (IQ-TREE)
6. Divergence time estimation (MCMCtree)

## Output

Results are saved in `results/timetree/`:

| File | Description |
|------|-------------|
| `supermatrix.phy` | Concatenated alignment |
| `species_tree.rooted.nwk` | Rooted ML tree |
| `timetree.final.nwk` | Final time tree (Newick) |
| `timetree.final.nex` | Final time tree (NEXUS, for FigTree) |

## Visualization

Open `results/timetree/timetree.final.nex` in FigTree:
1. Click "Node Labels" in left panel
2. Set "Display" to "date" to show estimated ages

## Expected Results

- ~1000+ single-copy orthologs (varies with species divergence)
- Land plants crown age: ~470 Mya
- Core angiosperms divergence: ~160 Mya
- Monocot-Dicot divergence: ~150 Mya
- Brassicaceae crown age: ~30 Mya
- Poaceae crown age: ~50 Mya
