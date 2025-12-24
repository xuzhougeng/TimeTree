"""
Rule: Run MCMCtree and post-process output (IQ2MC Step 3)
"""

rule run_mcmctree:
    """Run MCMCtree using IQ-TREE generated control file"""
    input:
        ctl=f"{config['output_dir']}/iq2mc/species.mcmctree.ctl",
        hessian=f"{config['output_dir']}/iq2mc/species.mcmctree.hessian"
    output:
        mcmc=f"{config['output_dir']}/mcmctree/mcmc.txt",
        figtree=f"{config['output_dir']}/mcmctree/FigTree.tre"
    params:
        ctl=f"{config['output_dir']}/iq2mc/species.mcmctree.ctl",
        outdir=f"{config['output_dir']}/mcmctree"
    log:
        f"{config['work_dir']}/logs/mcmctree.log"
    conda:
        "../envs/mcmctree.yaml"
    shell:
        """
        mkdir -p {params.outdir}
        cd {params.outdir}
        
        # Run mcmctree (requires modified PAML from iqtree/paml)
        mcmctree {params.ctl} 2>&1 | tee {log}
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

