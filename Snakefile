"""
TimeTree Snakemake Workflow
OrthoFinder SCO -> MSA -> Supermatrix -> IQ-TREE -> IQ2MC -> MCMCtree -> TimeTree
"""

from pathlib import Path

configfile: "config/config.yaml"

# Include rule modules
include: "workflow/rules/orthogroups.smk"
include: "workflow/rules/msa.smk"
include: "workflow/rules/supermatrix.smk"
include: "workflow/rules/iqtree.smk"
include: "workflow/rules/calibrate_tree.smk"
include: "workflow/rules/iq2mc.smk"
include: "workflow/rules/mcmctree.smk"
include: "workflow/rules/visualization.smk"

# Final target
rule all:
    input:
        f"{config['output_dir']}/timetree.final.nwk",
        f"{config['output_dir']}/timetree.final.nex",
        f"{config['output_dir']}/timetree.pdf",
        f"{config['output_dir']}/timetree.png",
        f"{config['output_dir']}/timetree.svg"

