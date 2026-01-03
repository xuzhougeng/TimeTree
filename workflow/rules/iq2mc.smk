"""
Rule: IQ-TREE --dating mcmctree to generate Hessian/ctl files (IQ2MC Step 2)
Supports multiple clock models (IND, CORR, EQUAL) running in parallel.
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
        hessian=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.mcmctree.hessian",
        ctl=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.mcmctree.ctl",
        rooted_nwk=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.rooted.nwk",
        dummy_aln=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.dummy.phy"
    params:
        prefix=lambda w: f"{config['output_dir']}/iq2mc/{w.clock_model}/species",
        model=config["iqtree"]["model"],
        use_partition=config["iqtree"]["use_partition"],
        clock=lambda w: w.clock_model,
        burnin=config["mcmctree"]["burnin"],
        samplefreq=config["mcmctree"]["samplefreq"],
        nsample=config["mcmctree"]["nsample"],
        birth=config["mcmctree"]["birth_rate"],
        death=config["mcmctree"]["death_rate"],
        sampling=config["mcmctree"]["sampling_fraction"],
        bin=config.get("binaries", {}).get("iqtree_dating", "")
    log:
        f"{config['work_dir']}/logs/iqtree_dating.{{clock_model}}.log"
    conda:
        "../envs/iqtree.yaml"
    threads: 8
    wildcard_constraints:
        clock_model="IND|CORR|EQUAL"
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
