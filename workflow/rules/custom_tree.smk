"""
Rule: Validate and copy custom species tree for use in downstream analysis.

When tree_method="custom" is specified, this rule:
1. Reads the user-provided Newick tree file
2. Validates that the tree is properly rooted
3. Validates that leaf names match species in the supermatrix
4. Copies the tree to the standard output location
"""


rule validate_custom_tree:
    """Validate and copy custom species tree"""
    input:
        tree=lambda wildcards: config.get("custom_tree", {}).get("tree_file", ""),
        supermatrix=f"{config['output_dir']}/supermatrix.faa"  # For species list validation
    output:
        rooted=f"{config['output_dir']}/species_tree.custom.rooted.nwk"
    log:
        f"{config['work_dir']}/logs/validate_custom_tree.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/process_custom_tree.py"
