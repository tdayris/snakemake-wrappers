"""
RScript provided
"""

rule ensembl_vep_haplotype_caller_ctc_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}_{nb}.vcf"],
    output:
        tsv=temp("vep/hc/{sample}_{version}_{manip}_{nb}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/hc/{sample}_{version}_{manip}_{nb}.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_hc.R")


rule ensembl_vep_haplotype_caller_wbc_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}_{version}_{manip}.vcf"],
    output:
        tsv=temp("vep/hc/{sample}_{version}_{manip}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/hc/{sample}_{version}_{manip}.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_hc.R")


rule ensembl_vep_haplotype_caller_baseline_postprocess:
    input:
        cancer_genes=config.get("cancer_genes", "Cancer.genes.cleaned.txt"),
        vcfs=["vep/annotate/{sample}.vcf"],
    output:
        tsv=temp("vep/hc/{sample}.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir=tmp,
    log:
        "logs/vep/hc/{sample}.log",
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
        str(workflow_source_dir / "scripts" / "ensemblVEP_hc.R")