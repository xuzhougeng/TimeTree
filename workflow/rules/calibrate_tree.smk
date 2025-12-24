"""
Rule: Annotate tree with fossil/node calibrations for MCMCtree
"""

rule annotate_calibrations:
    """Inject calibration constraints into rooted tree for MCMCtree"""
    input:
        tree=f"{config['output_dir']}/species_tree.rooted.nwk",
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

