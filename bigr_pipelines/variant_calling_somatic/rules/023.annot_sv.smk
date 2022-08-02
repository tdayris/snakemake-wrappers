rule annot_sv:
    input:
        vcf="facets/{sample}/{sample}.vcf.gz",
        vcf_idx="facets/{sample}/{sample}.vcf.gz.tbi",
    output:
        vcf=temp("annot_sv/{sample}.tsv"),
    threads: 1
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/annot_sv/raw/{sample}.log",
    params:
        install_dir=config["annot_sv"].get(
            "install_dir",
            "/mnt/beegfs/pipelines/snakemake-wrappers/bigr_pipelines/annot_sv/source/AnnotSV/",
        ),
        extra=config["annot_sv"].get(
            "extra", "-annotationMode both -snvIndelPASS 1 -tx ENSEMBL -SVinputInfo 0"
        ),
        genome=f"-genomeBuild {config['reference'].get('ncbi_build', 'GRCh38')}",
    conda:
        str(workflow_source_dir / "envs" / "annot_sv.yaml")
    shell:
        """
        export ANNOTSV="{params.install_dir}" &&
        "{params.install_dir}/bin/AnnotSV" {params.genome} -SVinputFile {input.vcf} -outputFile {output.vcf} > {log} 2>&1
        """


rule add_census:
    input:
        annot_sv="annot_sv/{sample}.tsv",
        census=config["reference"]["cancer_census"],
    output:
        census=protected("data_output/CNV/{sample}.tsv"),
    message:
        "Adding census data to {wildcards.sample}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 3,
        time_min=lambda wildcards, attempt: attempt * 5,
        tmpdir="tmp",
    log:
        "logs/annot_sv/census/{sample}.log",
    script:
        str(workflow_source_dir / "scripts" / "add_census.py")
