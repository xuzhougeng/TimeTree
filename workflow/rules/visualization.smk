"""
Rules: Visualize phylogenetic trees using ggtree
- Intermediate trees: rooted species tree, ASTRAL tree, IQ-TREE ML tree
- Final time tree with divergence times
"""


def get_intermediate_tree_plots(wildcards):
    """Return intermediate tree plots based on tree_method config."""
    method = config.get("tree_method", "concatenation")
    plots = []

    if method == "coalescent":
        # Coalescent method: ASTRAL trees
        plots.extend([
            f"{config['output_dir']}/plots/species_tree.astral.png",
            f"{config['output_dir']}/plots/species_tree.astral.rooted.png"
        ])
    elif method == "custom":
        # Custom method: only the user-provided rooted tree
        plots.extend([
            f"{config['output_dir']}/plots/species_tree.custom.rooted.png"
        ])
    else:
        # Concatenation method: IQ-TREE trees
        plots.extend([
            f"{config['output_dir']}/plots/species_tree.iqtree.png",
            f"{config['output_dir']}/plots/species_tree.rooted.png"
        ])

    return plots


# ============================================================================
# Aggregate Rules
# ============================================================================

rule plot_intermediate_trees:
    """Generate all intermediate tree visualizations"""
    input:
        get_intermediate_tree_plots
    output:
        touch(f"{config['output_dir']}/plots/.intermediate_done")


# ============================================================================
# Intermediate Tree Visualization
# ============================================================================

rule plot_rooted_species_tree:
    """Visualize the rooted species tree (concatenation method)"""
    input:
        tree=f"{config['output_dir']}/species_tree.rooted.nwk"
    output:
        pdf=f"{config['output_dir']}/plots/species_tree.rooted.pdf",
        png=f"{config['output_dir']}/plots/species_tree.rooted.png",
        svg=f"{config['output_dir']}/plots/species_tree.rooted.svg"
    params:
        title="Rooted Species Tree (Concatenation)",
        show_support=True
    log:
        f"{config['work_dir']}/logs/plot_rooted_species_tree.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_tree.R"


rule plot_iqtree_ml:
    """Visualize the IQ-TREE ML tree (unrooted)"""
    input:
        tree=f"{config['output_dir']}/iqtree/species.treefile"
    output:
        pdf=f"{config['output_dir']}/plots/species_tree.iqtree.pdf",
        png=f"{config['output_dir']}/plots/species_tree.iqtree.png",
        svg=f"{config['output_dir']}/plots/species_tree.iqtree.svg"
    params:
        title="IQ-TREE ML Tree (with Bootstrap)",
        show_support=True
    log:
        f"{config['work_dir']}/logs/plot_iqtree_ml.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_tree.R"


rule plot_astral_tree:
    """Visualize the ASTRAL species tree (unrooted)"""
    input:
        tree=f"{config['output_dir']}/coalescent/species.astral.nwk"
    output:
        pdf=f"{config['output_dir']}/plots/species_tree.astral.pdf",
        png=f"{config['output_dir']}/plots/species_tree.astral.png",
        svg=f"{config['output_dir']}/plots/species_tree.astral.svg"
    params:
        title="ASTRAL Species Tree (Coalescent)",
        show_support=True
    log:
        f"{config['work_dir']}/logs/plot_astral_tree.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_tree.R"


rule plot_rooted_astral_tree:
    """Visualize the rooted ASTRAL species tree"""
    input:
        tree=f"{config['output_dir']}/species_tree.astral.rooted.nwk"
    output:
        pdf=f"{config['output_dir']}/plots/species_tree.astral.rooted.pdf",
        png=f"{config['output_dir']}/plots/species_tree.astral.rooted.png",
        svg=f"{config['output_dir']}/plots/species_tree.astral.rooted.svg"
    params:
        title="Rooted ASTRAL Species Tree",
        show_support=True
    log:
        f"{config['work_dir']}/logs/plot_rooted_astral_tree.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_tree.R"


rule plot_custom_tree:
    """Visualize the custom rooted species tree"""
    input:
        tree=f"{config['output_dir']}/species_tree.custom.rooted.nwk"
    output:
        pdf=f"{config['output_dir']}/plots/species_tree.custom.rooted.pdf",
        png=f"{config['output_dir']}/plots/species_tree.custom.rooted.png",
        svg=f"{config['output_dir']}/plots/species_tree.custom.rooted.svg"
    params:
        title="Custom Species Tree (User-Provided)",
        show_support=False
    log:
        f"{config['work_dir']}/logs/plot_custom_tree.log"
    conda:
        "../envs/ggtree.yaml"
    script:
        "../scripts/plot_tree.R"


# ============================================================================
# Time Tree Visualization
# ============================================================================

rule plot_timetree:
    """Generate publication-quality time tree visualization with ggtree"""
    input:
        figtree=f"{config['output_dir']}/mcmctree/{{clock_model}}/FigTree.tre"
    output:
        pdf=f"{config['output_dir']}/timetree.{{clock_model}}.pdf",
        png=f"{config['output_dir']}/timetree.{{clock_model}}.png",
        svg=f"{config['output_dir']}/timetree.{{clock_model}}.svg"
    params:
        time_unit=config.get("visualization", {}).get("time_unit", "Ma"),
        clock_model=lambda w: w.clock_model
    log:
        f"{config['work_dir']}/logs/plot_timetree.{{clock_model}}.log"
    conda:
        "../envs/ggtree.yaml"
    wildcard_constraints:
        clock_model="IND|CORR|EQUAL"
    script:
        "../scripts/plot_timetree.R"
