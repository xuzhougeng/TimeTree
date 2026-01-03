# TimeTree Example: 8 Animal Species

This example demonstrates the TimeTree workflow using 8 animal species from Ensembl, representing major animal lineages from invertebrates to mammals.

## Species

| Species | Common Name | Lineage |
|---------|-------------|---------|
| Caenorhabditis elegans | Roundworm | Nematoda (outgroup) |
| Drosophila melanogaster | Fruit fly | Arthropoda (Insecta) |
| Danio rerio | Zebrafish | Actinopterygii (ray-finned fish) |
| Gallus gallus | Chicken | Aves (birds) |
| Bos taurus | Cow | Mammalia (Artiodactyla) |
| Canis lupus familiaris | Dog | Mammalia (Carnivora) |
| Mus musculus | Mouse | Mammalia (Rodentia) |
| Homo sapiens | Human | Mammalia (Primates) |

```newick
((Caenorhabditis_elegans,Drosophila_melanogaster),(Danio_rerio,(Gallus_gallus,((Bos_taurus,Canis_lupus_familiaris),(Mus_musculus,Homo_sapiens)))));
```

```tree
                    _____________________________ Caenorhabditis elegans
 __________________|
|                  |_____________________________ Drosophila melanogaster
|
|        _______________________________________ Danio rerio
|_______|
        |       ________________________________ Gallus gallus
        |______|
               |                 _____________________ Bos taurus
               |        ________|
               |       |        |_____________________ Canis lupus familiaris
               |_______|
                       |          _____________________ Mus musculus        
                       |________ |
                                 |_____________________ Homo sapiens

```

## Fossil Calibrations

| Calibration Point | Taxa | Age (Mya) | Reference |
|-------------------|------|-----------|-----------|
| Bilateria Crown | C. elegans, Human | 558-642 | Peterson et al. (2004) |
| Protostome-Deuterostome | Fruit fly, Human | 518-581 | Benton et al. (2015) |
| Osteichthyes Crown | Zebrafish, Human | 416-432 | Benton et al. (2015) |
| Amniota Crown | Chicken, Human | 312-330 | Benton et al. (2015) |
| Euarchontoglires Crown | Human, Mouse | 85-95 | Benton et al. (2015) |

## Directory Structure

```
example_animal/
├── calibrations.tsv          # Fossil calibration data
├── config.yaml               # Workflow configuration
├── data/
│   ├── cds/                  # CDS sequences
│   └── pep/                  # Protein sequences
├── orthofinder_input/        # OrthoFinder input (prepared)
├── orthofinder_analysis/     # OrthoFinder results
├── results/                  # TimeTree output
└── scripts/
    ├── download_animal_data.sh   # Download sequences from Ensembl
    ├── remove_isoforms.py        # Remove transcript isoforms
    ├── run_orthofinder.sh
    ├── setup_env.sh
    └── run_timetree.sh
```

## Quick Start

### Step 1: Download Data

```bash
./scripts/download_animal_data.sh
```

Downloads CDS and protein sequences from Ensembl (Release 113 for vertebrates, Ensembl Metazoa Release 60 for invertebrates), then automatically removes isoforms.

**What happens:**
1. Downloads `.cds.all.fa.gz` and `.pep.all.fa.gz` files from Ensembl
2. Removes isoforms using `remove_isoforms.py` - keeps only the longest transcript for each gene
3. Outputs uncompressed `.fa` files for faster downstream processing
4. Backs up original compressed files to `data/cds_original/` and `data/pep_original/`

### Step 2: Run OrthoFinder

```bash
# conda create -n of -c conda-forge -c bioconda orthofinder
./scripts/run_orthofinder.sh
```

Identifies orthogroups and single-copy orthologs. Update `config.yaml` with the correct `results_dir` path after OrthoFinder completes.

### Step 3: Setup Environment

```bash
./scripts/setup_env.sh
```

Creates `timetree` conda environment with all required tools.

### Step 4: Run TimeTree

```bash
./scripts/run_timetree.sh
```

Executes the complete workflow:
1. Extract single-copy orthologs
2. Multiple sequence alignment (MAFFT)
3. Alignment trimming (trimAl)
4. Build gene trees (IQ-TREE)
5. Species tree inference (ASTRAL coalescent method)
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

- ~2000+ single-copy orthologs across all 8 species
- Bilateria crown age: ~600 Mya
- Vertebrate-Invertebrate divergence: ~550 Mya
- Fish-Tetrapod divergence: ~420 Mya
- Amniota (Bird-Mammal) divergence: ~320 Mya
- Primate-Rodent divergence: ~90 Mya

## Notes

- Vertebrate genomes are significantly larger than plant genomes; OrthoFinder may take longer
- The species selection covers key nodes in animal evolution for robust dating
- C. elegans serves as a distant outgroup for rooting the tree
