"""
Rule: Run MCMCtree and post-process output (IQ2MC Step 3)
Supports multiple clock models (IND, CORR, EQUAL) running in parallel.
"""

rule run_mcmctree:
    """Run MCMCtree using IQ-TREE generated control file"""
    input:
        ctl=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.mcmctree.ctl",
        hessian=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.mcmctree.hessian",
        rooted_nwk=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.rooted.nwk",
        dummy_aln=f"{config['output_dir']}/iq2mc/{{clock_model}}/species.dummy.phy"
    output:
        mcmc=f"{config['output_dir']}/mcmctree/{{clock_model}}/mcmc.txt",
        figtree=f"{config['output_dir']}/mcmctree/{{clock_model}}/FigTree.tre"
    params:
        outdir=lambda w: f"{config['output_dir']}/mcmctree/{w.clock_model}",
        mcmctree_bin=config.get("binaries", {}).get("mcmctree", "mcmctree"),
        parse_ctl_script=str(Path(workflow.basedir) / "workflow/scripts/parse_mcmctree_ctl.py"),
        localize_ctl_script=str(Path(workflow.basedir) / "workflow/scripts/localize_mcmctree_ctl.py"),
        clean_tree_script=str(Path(workflow.basedir) / "workflow/scripts/clean_mcmctree_tree.py"),
        adjust_rootage_script=str(Path(workflow.basedir) / "workflow/scripts/adjust_rootage.py"),
        ensure_print_script=str(Path(workflow.basedir) / "workflow/scripts/ensure_mcmctree_print.py")
    log:
        f"{config['work_dir']}/logs/mcmctree.{{clock_model}}.log"
    conda:
        "../envs/mcmctree.yaml"
    wildcard_constraints:
        clock_model="IND|CORR|EQUAL"
    shell:
        """
        mkdir -p {params.outdir}
        LOG_FILE="$(pwd)/{log}"
        OUTDIR_ABS="$(pwd)/{params.outdir}"
        mkdir -p "$(dirname "$LOG_FILE")"
        PARSE_CTL_SCRIPT="{params.parse_ctl_script}"
        LOCALIZE_CTL_SCRIPT="{params.localize_ctl_script}"
        CLEAN_TREE_SCRIPT="{params.clean_tree_script}"
        ADJUST_ROOTAGE_SCRIPT="{params.adjust_rootage_script}"
        ENSURE_PRINT_SCRIPT="{params.ensure_print_script}"

        # MCMCtree ctl uses relative paths (seqfile/treefile). Run inside a staging dir
        # containing the ctl + referenced files (as recommended in IQ-TREE Dating docs).
        RUNDIR="{params.outdir}/run"
        mkdir -p "$RUNDIR"

        cp -f {input.ctl} {input.hessian} {input.rooted_nwk} {input.dummy_aln} "$RUNDIR/"
        cd "$RUNDIR"

        # Run mcmctree (IQ2MC step3)
        python "$LOCALIZE_CTL_SCRIPT" species.mcmctree.ctl
        python "$CLEAN_TREE_SCRIPT" species.rooted.nwk
        python "$ADJUST_ROOTAGE_SCRIPT" species.mcmctree.ctl species.rooted.nwk
        python "$ENSURE_PRINT_SCRIPT" species.mcmctree.ctl
        MCMCFILE=$(python "$PARSE_CTL_SCRIPT" species.mcmctree.ctl)
        if [ -z "$MCMCFILE" ]; then
            MCMCFILE="mcmc.txt"
        fi

        {params.mcmctree_bin} species.mcmctree.ctl 2>&1 | tee "$LOG_FILE"

        # Collect key outputs
        if [ -f "$MCMCFILE" ]; then
            cp -f "$MCMCFILE" "$OUTDIR_ABS/mcmc.txt"
        else
            echo "WARNING: mcmcfile not found: $MCMCFILE" >&2
        fi
        cp -f FigTree.tre "$OUTDIR_ABS/FigTree.tre"
        """


rule postprocess_timetree:
    """Post-process MCMCtree output to final time tree formats"""
    input:
        figtree=f"{config['output_dir']}/mcmctree/{{clock_model}}/FigTree.tre",
        mcmc=f"{config['output_dir']}/mcmctree/{{clock_model}}/mcmc.txt"
    output:
        nwk=f"{config['output_dir']}/timetree.{{clock_model}}.final.nwk",
        nex=f"{config['output_dir']}/timetree.{{clock_model}}.final.nex"
    log:
        f"{config['work_dir']}/logs/postprocess_timetree.{{clock_model}}.log"
    conda:
        "../envs/py.yaml"
    wildcard_constraints:
        clock_model="IND|CORR|EQUAL"
    script:
        "../scripts/mcmctree_postprocess.py"
