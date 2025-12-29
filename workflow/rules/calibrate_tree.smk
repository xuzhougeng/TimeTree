"""
Rule: Annotate tree with fossil/node calibrations for MCMCtree
"""


def get_species_tree_for_calibration(wildcards):
    """Get the rooted species tree based on tree_method config."""
    method = config.get("tree_method", "concatenation")
    if method == "coalescent":
        return f"{config['output_dir']}/species_tree.astral.rooted.nwk"
    else:  # concatenation (default)
        return f"{config['output_dir']}/species_tree.rooted.nwk"


rule annotate_calibrations:
    """Inject calibration constraints into rooted tree for MCMCtree"""
    input:
        tree=get_species_tree_for_calibration,
        calibrations=config["calibration_table"]
    output:
        calibrated=f"{config['output_dir']}/species_tree.calibrated.nwk",
        log_file=f"{config['output_dir']}/calibration_mapping.log"
    params:
        unit=config["calibration_unit"]
    log:
        f"{config['work_dir']}/logs/annotate_calibrations.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/annotate_calibrations.py"

