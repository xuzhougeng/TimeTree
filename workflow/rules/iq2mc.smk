"""
Rule: IQ-TREE --dating mcmctree to generate Hessian/ctl files (IQ2MC Step 2)
"""


def get_iq2mc_inputs(wildcards):
    """Get inputs for iqtree_dating_mcmctree based on tree_method."""
    inputs = {
        "supermatrix": f"{config['output_dir']}/supermatrix.phy",
        "partitions": f"{config['output_dir']}/partitions.nex",
        "tree": f"{config['output_dir']}/species_tree.calibrated.nwk",
    }
    # Only require iqtree_done for concatenation method
    # (coalescent and custom methods skip IQ-TREE ML tree inference)
    if config.get("tree_method", "concatenation") == "concatenation":
        inputs["iqtree_done"] = f"{config['output_dir']}/iqtree/species.iqtree"
    return inputs


rule iqtree_dating_mcmctree:
    """Generate Hessian and control file for MCMCtree"""
    input:
        unpack(get_iq2mc_inputs)
    output:
        hessian=f"{config['output_dir']}/iq2mc/species.mcmctree.hessian",
        ctl=f"{config['output_dir']}/iq2mc/species.mcmctree.ctl",
        rooted_nwk=f"{config['output_dir']}/iq2mc/species.rooted.nwk",
        dummy_aln=f"{config['output_dir']}/iq2mc/species.dummy.phy"
    params:
        prefix=f"{config['output_dir']}/iq2mc/species",
        model=config["iqtree"]["model"],
        use_partition=config["iqtree"]["use_partition"],
        clock=config["mcmctree"]["clock_model"],
        burnin=config["mcmctree"]["burnin"],
        samplefreq=config["mcmctree"]["samplefreq"],
        nsample=config["mcmctree"]["nsample"],
        birth=config["mcmctree"]["birth_rate"],
        death=config["mcmctree"]["death_rate"],
        sampling=config["mcmctree"]["sampling_fraction"],
        bin=config.get("binaries", {}).get("iqtree_dating", "")
    log:
        f"{config['work_dir']}/logs/iqtree_dating.log"
    conda:
        "../envs/iqtree.yaml"
    threads: 8
    shell:
        """
        # IQ2MC step2 (Dating doc): REQUIRE iqtree3 for --dating mcmctree
        mkdir -p $(dirname {params.prefix})

        IQTREE_DATING_BIN="{params.bin}"
        if [ -z "$IQTREE_DATING_BIN" ]; then
            IQTREE_DATING_BIN="iqtree3"
        fi

        if ! command -v "$IQTREE_DATING_BIN" >/dev/null 2>&1; then
            echo "ERROR: IQ2MC step2 requires IQ-TREE 3 (iqtree3). Not found: $IQTREE_DATING_BIN" >&2
            echo "Set config.yaml: binaries.iqtree_dating to the absolute path of iqtree3, or add iqtree3 to PATH." >&2
            exit 1
        fi

        MODEL="{params.model}"
        
        PARTITION_OPT=""
        if [ "{params.use_partition}" = "True" ]; then
            PARTITION_OPT="-p {input.partitions}"
        fi
        
        $IQTREE_DATING_BIN -s {input.supermatrix} $PARTITION_OPT \
            -m "$MODEL" \
            -te {input.tree} \
            --dating mcmctree \
            --mcmc-iter {params.burnin},{params.samplefreq},{params.nsample} \
            --mcmc-bds {params.birth},{params.death},{params.sampling} \
            --mcmc-clock {params.clock} \
            -T {threads} \
            --prefix {params.prefix} \
            -redo \
            2>&1 | tee {log}
        """
