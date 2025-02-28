"""
R script provided
"""

rule postprocess_ensembl_vep_mutect_ctc:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}_{nb}.mutect.vcf"],
    output:
        tsv=temp("vep/mutect/{sample}_{version}_{manip}_{nb}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/mutect/{sample}_{version}_{manip}_{nb}.log",
    params:
        organism=config.get("vep_db", "hg38"),
    container:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / ".."
            / "singularity"
            / "mambaforge_4.14.0-0.sif"
        )
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        str(workflow_source_dir / "scripts" / "ensemblVEP_mutect.R")



rule postprocess_ensembl_vep_mutect_wbc:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}.mutect.vcf"],
    output:
        tsv=temp("vep/mutect/{sample}_{version}_{manip}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/mutect/{sample}_{version}_{manip}.log",
    params:
        organism=config.get("vep_db", "hg38"),
    container:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / ".."
            / "singularity"
            / "mambaforge_4.14.0-0.sif"
        )
    conda:
        str(workflow_source_dir / "envs" / "r.yaml")
    script:
        str(workflow_source_dir / "scripts" / "ensemblVEP_mutect.R")