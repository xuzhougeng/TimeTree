"""
Rule: Run MCMCtree and post-process output (IQ2MC Step 3)
"""

rule run_mcmctree:
    """Run MCMCtree using IQ-TREE generated control file"""
    input:
        ctl=f"{config['output_dir']}/iq2mc/species.mcmctree.ctl",
        hessian=f"{config['output_dir']}/iq2mc/species.mcmctree.hessian",
        rooted_nwk=f"{config['output_dir']}/iq2mc/species.rooted.nwk",
        dummy_aln=f"{config['output_dir']}/iq2mc/species.dummy.phy"
    output:
        mcmc=f"{config['output_dir']}/mcmctree/mcmc.txt",
        figtree=f"{config['output_dir']}/mcmctree/FigTree.tre"
    params:
        outdir=f"{config['output_dir']}/mcmctree",
        mcmctree_bin=config.get("binaries", {}).get("mcmctree", "mcmctree")
    log:
        f"{config['work_dir']}/logs/mcmctree.log"
    conda:
        "../envs/mcmctree.yaml"
    shell:
        """
        mkdir -p {params.outdir}

        # MCMCtree ctl uses relative paths (seqfile/treefile). Run inside a staging dir
        # containing the ctl + referenced files (as recommended in IQ-TREE Dating docs).
        RUNDIR="{params.outdir}/run"
        mkdir -p "$RUNDIR"

        cp -f {input.ctl} {input.hessian} {input.rooted_nwk} {input.dummy_aln} "$RUNDIR/"
        cd "$RUNDIR"

        # Run mcmctree (IQ2MC step3)
        {params.mcmctree_bin} species.mcmctree.ctl 2>&1 | tee {log}

        # Collect key outputs
        cp -f mcmc.txt "{params.outdir}/mcmc.txt"
        cp -f FigTree.tre "{params.outdir}/FigTree.tre"
        """


rule postprocess_timetree:
    """Post-process MCMCtree output to final time tree formats"""
    input:
        figtree=f"{config['output_dir']}/mcmctree/FigTree.tre",
        mcmc=f"{config['output_dir']}/mcmctree/mcmc.txt"
    output:
        nwk=f"{config['output_dir']}/timetree.final.nwk",
        nex=f"{config['output_dir']}/timetree.final.nex"
    log:
        f"{config['work_dir']}/logs/postprocess_timetree.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/mcmctree_postprocess.py"

