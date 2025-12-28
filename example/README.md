# TimeTree Example: 5 Plant Species

This example demonstrates the TimeTree workflow using 5 plant species from Ensembl Plants.

## Species

| Species | Common Name |
|---------|-------------|
| Arabidopsis thaliana | Thale cress |
| Oryza sativa | Rice |
| Zea mays | Maize |
| Solanum lycopersicum | Tomato |
| Vitis vinifera | Grape |

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
    ├── run_orthofinder.sh
    ├── setup_env.sh
    └── run_timetree.sh
```

## Quick Start

### Step 1: Download Data

```bash
./scripts/download_plant_data.sh
```

Downloads CDS and protein sequences from Ensembl Plants Release 60.

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

- ~460 single-copy orthologs
- Monocot-Dicot divergence: ~150 Mya
- Poaceae crown age: ~50 Mya
