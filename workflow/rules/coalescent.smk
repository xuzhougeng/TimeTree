"""
Rule: Coalescent-based species tree inference using gene trees and ASTRAL
This method is more robust to Incomplete Lineage Sorting (ILS) than concatenation.

Workflow:
1. Build individual gene trees for each orthogroup (IQ-TREE)
2. Concatenate gene trees into a single file
3. Run ASTRAL to infer species tree from gene trees
4. Root the species tree
"""

from pathlib import Path


rule build_gene_tree:
    """Build ML gene tree for each orthogroup using IQ-TREE"""
    input:
        aln=f"{config['work_dir']}/msa/{{og}}.aln.faa"
    output:
        treefile=f"{config['work_dir']}/gene_trees/{{og}}.treefile"
    params:
        model=config.get("coalescent", {}).get("gene_tree_model", "MFP"),
        prefix=f"{config['work_dir']}/gene_trees/{{og}}",
        bin=config.get("binaries", {}).get("iqtree_ml", "")
    log:
        f"{config['work_dir']}/logs/gene_trees/{{og}}.log"
    conda:
        "../envs/iqtree.yaml"
    threads: config.get("coalescent", {}).get("gene_tree_threads", 2)
    shell:
        """
        mkdir -p $(dirname {params.prefix})

        IQTREE_BIN="{params.bin}"
        if [ -z "$IQTREE_BIN" ]; then
            if command -v iqtree3 >/dev/null 2>&1; then
                IQTREE_BIN="iqtree3"
            elif command -v iqtree >/dev/null 2>&1; then
                IQTREE_BIN="iqtree"
            else
                echo "ERROR: cannot find IQ-TREE binary (iqtree3/iqtree) in PATH" >&2
                exit 1
            fi
        fi

        $IQTREE_BIN -s {input.aln} \
            -m {params.model} \
            -T {threads} \
            --prefix {params.prefix} \
            -redo \
            2>&1 | tee {log}
        """


def get_gene_trees(wildcards):
    """Get all gene tree files based on SCO checkpoint"""
    checkpoint_output = checkpoints.extract_sco.get(**wildcards).output[0]
    sco_dir = Path(checkpoint_output)
    ogs = [f.stem for f in sco_dir.glob("OG*.faa")]
    return expand(f"{config['work_dir']}/gene_trees/{{og}}.treefile", og=ogs)


rule concat_gene_trees:
    """Concatenate all gene trees into a single file for ASTRAL input"""
    input:
        trees=get_gene_trees
    output:
        combined=f"{config['output_dir']}/coalescent/gene_trees.nwk"
    log:
        f"{config['work_dir']}/logs/concat_gene_trees.log"
    run:
        import os
        os.makedirs(os.path.dirname(output.combined), exist_ok=True)

        with open(output.combined, 'w') as out:
            for tree_file in input.trees:
                with open(tree_file) as f:
                    tree = f.read().strip()
                    if tree:
                        out.write(tree + '\n')

        with open(log[0], 'w') as logf:
            logf.write(f"Combined {len(input.trees)} gene trees into {output.combined}\n")


rule run_astral:
    """Run ASTRAL to infer consensus species tree from gene trees

    ASTRAL uses the multi-species coalescent model to infer a consensus
    species tree from a collection of gene trees. It is statistically
    consistent under the coalescent model, making it robust to ILS.

    The output is a consensus species tree that minimizes quartet distance
    across all input gene trees. Use astral_opts in config to add options:
    - "-t 2" to output quartet support values as branch labels
    - "-t 8" to output local posterior probabilities
    """
    input:
        gene_trees=f"{config['output_dir']}/coalescent/gene_trees.nwk"
    output:
        species_tree=f"{config['output_dir']}/coalescent/species.astral.nwk",
        log_file=f"{config['output_dir']}/coalescent/astral.log"
    params:
        extra_opts=config.get("coalescent", {}).get("astral_opts", ""),
        bin=config.get("coalescent", {}).get("astral_bin", "astral")
    log:
        f"{config['work_dir']}/logs/astral.log"
    conda:
        "../envs/phylo.yaml"
    threads: config.get("coalescent", {}).get("astral_threads", 4)
    shell:
        """
        mkdir -p $(dirname {output.species_tree})

        # Check for ASTRAL binary (astral or astral.jar)
        ASTRAL_CMD=""
        if [ -n "{params.bin}" ] && [ "{params.bin}" != "astral" ]; then
            # User specified a custom path
            if [[ "{params.bin}" == *.jar ]]; then
                ASTRAL_CMD="java -jar {params.bin}"
            else
                ASTRAL_CMD="{params.bin}"
            fi
        else
            # Auto-detect ASTRAL
            if command -v astral >/dev/null 2>&1; then
                ASTRAL_CMD="astral"
            elif command -v astral-pro >/dev/null 2>&1; then
                ASTRAL_CMD="astral-pro"
            elif [ -f "$CONDA_PREFIX/share/astral-tree/astral.jar" ]; then
                ASTRAL_CMD="java -jar $CONDA_PREFIX/share/astral-tree/astral.jar"
            else
                echo "ERROR: Cannot find ASTRAL binary. Install with 'conda install -c bioconda astral-tree' or specify coalescent.astral_bin in config." >&2
                exit 1
            fi
        fi

        echo "Using ASTRAL command: $ASTRAL_CMD" | tee {log}

        $ASTRAL_CMD -i {input.gene_trees} -o {output.species_tree} \
            {params.extra_opts} \
            2>&1 | tee -a {log}

        cp {log} {output.log_file}
        """


rule root_astral_tree:
    """Root the ASTRAL species tree using outgroup (required for coalescent method)"""
    input:
        tree=f"{config['output_dir']}/coalescent/species.astral.nwk"
    output:
        rooted=f"{config['output_dir']}/species_tree.astral.rooted.nwk"
    params:
        outgroup=config["outgroup_taxa"],
        require_outgroup=True
    log:
        f"{config['work_dir']}/logs/root_astral_tree.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/root_tree.py"


# Summary statistics for gene trees
rule gene_tree_stats:
    """Generate summary statistics for gene trees"""
    input:
        gene_trees=f"{config['output_dir']}/coalescent/gene_trees.nwk"
    output:
        stats=f"{config['output_dir']}/coalescent/gene_tree_stats.txt"
    log:
        f"{config['work_dir']}/logs/gene_tree_stats.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/gene_tree_stats.py"
