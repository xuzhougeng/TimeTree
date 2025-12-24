"""
Rule: Multiple sequence alignment and trimming for each orthogroup
"""

rule mafft_align:
    """Align sequences with MAFFT"""
    input:
        fasta=f"{config['work_dir']}/orthogroups_sco/{{og}}.faa"
    output:
        aln=f"{config['work_dir']}/msa/{{og}}.raw.aln"
    params:
        opts=config["msa"]["mafft_opts"]
    log:
        f"{config['work_dir']}/logs/mafft/{{og}}.log"
    conda:
        "../envs/phylo.yaml"
    threads: 2
    shell:
        "mafft {params.opts} --thread {threads} {input.fasta} > {output.aln} 2> {log}"


rule trim_alignment:
    """Trim alignment with trimal or clipkit"""
    input:
        aln=f"{config['work_dir']}/msa/{{og}}.raw.aln"
    output:
        trimmed=f"{config['work_dir']}/msa/{{og}}.aln.faa"
    params:
        method=config["trim"]["method"],
        trimal_opts=config["trim"]["trimal_opts"]
    log:
        f"{config['work_dir']}/logs/trim/{{og}}.log"
    conda:
        "../envs/phylo.yaml"
    shell:
        """
        if [ "{params.method}" = "trimal" ]; then
            trimal -in {input.aln} -out {output.trimmed} {params.trimal_opts} 2> {log}
        else
            clipkit {input.aln} -o {output.trimmed} 2> {log}
        fi
        """


def get_trimmed_alignments(wildcards):
    """Get all trimmed alignment files"""
    checkpoint_output = checkpoints.extract_sco.get(**wildcards).output[0]
    sco_dir = Path(checkpoint_output)
    ogs = [f.stem for f in sco_dir.glob("OG*.faa")]
    return expand(f"{config['work_dir']}/msa/{{og}}.aln.faa", og=ogs)

