"""
Rscript provided
"""

rule ensembl_vep_bcr_ctc_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}_{nb}.brc.vcf"],
    output:
        tsv=temp("vep/bcr/{sample}_{version}_{manip}_{nb}.tsv"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/bcr/{sample}_{version}_{manip}_{nb}.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_bcr.R")



rule ensembl_vep_bcr_ctc_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}.brc.vcf"],
    output:
        tsv=temp("vep/bcr/{sample}_{version}_{manip}.tsv"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/bcr/{sample}_{version}_{manip}.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_bcr.R")


rule ensembl_vep_bcr_baseline_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}.baseline.brc.vcf"],
    output:
        tsv=temp("vep/bcr/{sample}.baseline.tsv"),
    threads: 1
    resources:
        mem_mb=get_20gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/bcr/{sample}.baseline.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_bcr.R")