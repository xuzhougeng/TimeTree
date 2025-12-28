"""
Rule: Visualize time-calibrated tree using ggtree
"""

rule plot_timetree:
    """Generate publication-quality time tree visualization with ggtree"""
    input:
        figtree=f"{config['output_dir']}/mcmctree/FigTree.tre"
    output:
        pdf=f"{config['output_dir']}/timetree.pdf",
        png=f"{config['output_dir']}/timetree.png",
        svg=f"{config['output_dir']}/timetree.svg"
    params:
        time_unit=config.get("visualization", {}).get("time_unit", "Ma")
    log:
        f"{config['work_dir']}/logs/plot_timetree.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_timetree.R"
