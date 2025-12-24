"""
Rule: IQ-TREE --dating mcmctree to generate Hessian/ctl files (IQ2MC Step 2)
"""

rule iqtree_dating_mcmctree:
    """Generate Hessian and control file for MCMCtree"""
    input:
        supermatrix=f"{config['output_dir']}/supermatrix.phy",
        partitions=f"{config['output_dir']}/partitions.nex",
        tree=f"{config['output_dir']}/species_tree.calibrated.nwk",
        iqtree_done=f"{config['output_dir']}/iqtree/species.iqtree"  # ensure model is fitted
    output:
        hessian=f"{config['output_dir']}/iq2mc/species.mcmctree.hessian",
        ctl=f"{config['output_dir']}/iq2mc/species.mcmctree.ctl",
        rooted_nwk=f"{config['output_dir']}/iq2mc/species.rooted.nwk",
        dummy_aln=f"{config['output_dir']}/iq2mc/species.dummy.phy"
    params:
        prefix=f"{config['output_dir']}/iq2mc/species",
        use_partition=config["iqtree"]["use_partition"],
        clock=config["mcmctree"]["clock_model"],
        burnin=config["mcmctree"]["burnin"],
        samplefreq=config["mcmctree"]["samplefreq"],
        nsample=config["mcmctree"]["nsample"],
        birth=config["mcmctree"]["birth_rate"],
        death=config["mcmctree"]["death_rate"],
        sampling=config["mcmctree"]["sampling_fraction"]
    log:
        f"{config['work_dir']}/logs/iqtree_dating.log"
    conda:
        "../envs/iqtree.yaml"
    threads: 8
    shell:
        """
        # Get model from previous IQ-TREE run
        MODEL=$(grep "Best-fit model" {input.iqtree_done} | head -1 | awk '{{print $NF}}' || echo "LG+G4")
        
        PARTITION_OPT=""
        if [ "{params.use_partition}" = "True" ]; then
            PARTITION_OPT="-p {input.partitions}"
        fi
        
        iqtree2 -s {input.supermatrix} $PARTITION_OPT \
            -m $MODEL \
            -te {input.tree} \
            --dating mcmctree \
            --mcmc-iter {params.burnin},{params.samplefreq},{params.nsample} \
            --mcmc-bds {params.birth},{params.death},{params.sampling} \
            --mcmc-clock {params.clock} \
            -T {threads} \
            --prefix {params.prefix} \
            2>&1 | tee {log}
        """

