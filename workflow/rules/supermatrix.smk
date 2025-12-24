"""
Rule: Concatenate alignments into supermatrix with partition file
"""

rule concat_supermatrix:
    """Concatenate all trimmed alignments into supermatrix"""
    input:
        alignments=get_trimmed_alignments
    output:
        supermatrix=f"{config['output_dir']}/supermatrix.phy",
        partitions=f"{config['output_dir']}/partitions.nex",
        fasta=f"{config['output_dir']}/supermatrix.faa"
    log:
        f"{config['work_dir']}/logs/concat_supermatrix.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/concat_alignments.py"

