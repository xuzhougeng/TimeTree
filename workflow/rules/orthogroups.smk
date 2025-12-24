"""
Rule: Extract single-copy orthogroups from OrthoFinder results
"""

checkpoint extract_sco:
    """Extract SCO FASTA files from OrthoFinder output"""
    input:
        sco_list=f"{config['results_dir']}/Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
        og_dir=f"{config['results_dir']}/Orthogroup_Sequences",
        proteomes=config["proteomes_dir"]
    output:
        directory(f"{config['work_dir']}/orthogroups_sco"),
        log_file=f"{config['work_dir']}/orthogroups_sco/extraction.log"
    params:
        min_taxa=config["sco"]["min_taxa"],
        allow_missing=config["sco"]["allow_missing"]
    log:
        f"{config['work_dir']}/logs/extract_sco.log"
    conda:
        "../envs/py.yaml"
    script:
        "../scripts/extract_sco_fastas.py"


def get_sco_fastas(wildcards):
    """Get list of SCO FASTA files after checkpoint"""
    checkpoint_output = checkpoints.extract_sco.get(**wildcards).output[0]
    sco_dir = Path(checkpoint_output)
    return sorted(sco_dir.glob("OG*.faa"))

