"""
TimeTree Snakemake Workflow
OrthoFinder SCO -> MSA -> Species Tree -> IQ2MC -> MCMCtree -> TimeTree

Species tree methods:
  - concatenation: Supermatrix + IQ-TREE ML tree (traditional)
                   Outgroup optional (midpoint rooting if not specified)
  - coalescent: Gene trees + ASTRAL consensus tree (robust to ILS)
                REQUIRES outgroup setting for rooting
  - custom: User-provided pre-built rooted tree (skip tree inference)
            REQUIRES custom_tree.tree_file to be set
"""

from pathlib import Path

configfile: "config/config.yaml"

# ============================================================================
# Clock Models Configuration
# ============================================================================

def get_clock_models():
    """Get clock models to run with backward compatibility."""
    mcmctree_config = config.get("mcmctree", {})
    if "clock_models" in mcmctree_config:
        return mcmctree_config["clock_models"]
    elif "clock_model" in mcmctree_config:
        # Backward compatibility: single model as list
        return [mcmctree_config["clock_model"]]
    else:
        return ["IND"]  # Default

CLOCK_MODELS = get_clock_models()

# ============================================================================
# Configuration Validation
# ============================================================================

def validate_config():
    """Validate configuration settings."""
    tree_method = config.get("tree_method", "concatenation")
    valid_methods = ["concatenation", "coalescent", "custom"]

    if tree_method not in valid_methods:
        raise ValueError(
            f"Invalid tree_method: '{tree_method}'. "
            f"Must be one of: {valid_methods}"
        )

    # Validate coalescent method requirements
    if tree_method == "coalescent":
        outgroup = config.get("outgroup_taxa", [])
        if not outgroup:
            raise ValueError(
                "tree_method='coalescent' requires outgroup_taxa to be set.\n"
                "Please specify at least one outgroup species in config.yaml:\n"
                "  outgroup_taxa: ['species_name']"
            )

    # Validate custom method requirements
    if tree_method == "custom":
        custom_tree = config.get("custom_tree", {})
        tree_file = custom_tree.get("tree_file", "")
        if not tree_file:
            raise ValueError(
                "tree_method='custom' requires custom_tree.tree_file to be set.\n"
                "Please specify the path to your rooted species tree in config.yaml:\n"
                "  custom_tree:\n"
                "    tree_file: 'path/to/your/species_tree.nwk'"
            )
        # Check if file exists
        if not Path(tree_file).exists():
            raise ValueError(
                f"Custom tree file not found: {tree_file}\n"
                "Please check the path in custom_tree.tree_file"
            )

# Run validation at workflow start
validate_config()

# Include rule modules
include: "workflow/rules/orthogroups.smk"
include: "workflow/rules/msa.smk"
include: "workflow/rules/supermatrix.smk"
include: "workflow/rules/iqtree.smk"
include: "workflow/rules/coalescent.smk"
include: "workflow/rules/custom_tree.smk"
include: "workflow/rules/calibrate_tree.smk"
include: "workflow/rules/iq2mc.smk"
include: "workflow/rules/mcmctree.smk"
include: "workflow/rules/visualization.smk"


def get_rooted_species_tree(wildcards):
    """Return the rooted species tree based on tree_method config."""
    method = config.get("tree_method", "concatenation")
    if method == "coalescent":
        return f"{config['output_dir']}/species_tree.astral.rooted.nwk"
    elif method == "custom":
        return f"{config['output_dir']}/species_tree.custom.rooted.nwk"
    else:  # concatenation (default)
        return f"{config['output_dir']}/species_tree.rooted.nwk"


# Final target
rule all:
    input:
        expand(f"{config['output_dir']}/timetree.{{clock}}.final.nwk", clock=CLOCK_MODELS),
        expand(f"{config['output_dir']}/timetree.{{clock}}.final.nex", clock=CLOCK_MODELS),
        expand(f"{config['output_dir']}/timetree.{{clock}}.pdf", clock=CLOCK_MODELS),
        expand(f"{config['output_dir']}/timetree.{{clock}}.png", clock=CLOCK_MODELS),
        expand(f"{config['output_dir']}/timetree.{{clock}}.svg", clock=CLOCK_MODELS)

