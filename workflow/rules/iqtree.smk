"""
Rule: IQ-TREE ML tree inference (IQ2MC Step 1)
"""

rule iqtree_ml:
    """Infer ML tree with IQ-TREE"""
    input:
        supermatrix=f"{config['output_dir']}/supermatrix.phy",
        partitions=f"{config['output_dir']}/partitions.nex"
    output:
        treefile=f"{config['output_dir']}/iqtree/species.treefile",
        iqtree=f"{config['output_dir']}/iqtree/species.iqtree",
        log_file=f"{config['output_dir']}/iqtree/species.log"
    params:
        model=config["iqtree"]["model"],
        use_partition=config["iqtree"]["use_partition"],
        bootstrap=config["iqtree"]["bootstrap"],
        threads=config["iqtree"]["threads"],
        outgroup=config["outgroup_taxa"],
        prefix=f"{config['output_dir']}/iqtree/species",
        bin=config.get("binaries", {}).get("iqtree_ml", "")
    log:
        f"{config['work_dir']}/logs/iqtree_ml.log"
    conda:
        "../envs/iqtree.yaml"
    threads: config["iqtree"].get("threads", 8) if isinstance(config["iqtree"].get("threads"), int) else 8
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

        PARTITION_OPT=""
        if [ "{params.use_partition}" = "True" ]; then
            PARTITION_OPT="-p {input.partitions}"
        fi
        
        BOOT_OPT=""
        if [ "{params.bootstrap}" -gt 0 ]; then
            BOOT_OPT="-B {params.bootstrap}"
        fi
        
        OUTGROUP_OPT=""
        if [ -n "{params.outgroup}" ] && [ "{params.outgroup}" != "[]" ]; then
            # Convert list to comma-separated string
            OG=$(echo "{params.outgroup}" | tr -d "[]'" | tr ' ' ',')
            if [ -n "$OG" ]; then
                OUTGROUP_OPT="-o $OG"
            fi
        fi
        
        $IQTREE_BIN -s {input.supermatrix} $PARTITION_OPT \
            -m {params.model} $BOOT_OPT $OUTGROUP_OPT \
            -T {params.threads} \
            --prefix {params.prefix} \
            2>&1 | tee {log}
        """


rule root_tree:
    """Ensure tree is rooted (midpoint if no outgroup specified)"""
    input:
        tree=f"{config['output_dir']}/iqtree/species.treefile"
    output:
        rooted=f"{config['output_dir']}/species_tree.rooted.nwk"
    params:
        outgroup=config["outgroup_taxa"]
    log:
        f"{config['work_dir']}/logs/root_tree.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/root_tree.py"
