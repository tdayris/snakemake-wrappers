###############################
### Prepare nf-core chipseq ###
###############################

# Build nf-core design file
rule nextflox_design:
    input:
        fastq_files=expand(
            "reads/{sample}.{stream}.fq.gz",
            sample=design.Sample_id,
            stream=["1", "2"]
        )
    output:
        "design.csv"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        time_min=lambda wildcards, attempt: attempt * 5,
        tmpdir="tmp"
    log:
        "logs/nf-core/design.log"
    env:
        "envs/python.yaml"
    params:
        design=design.copy(),
        fastq_links=fastq_links.copy()
    script:
        "scripts/001.design_to_nfcore.py"


# Build nf-core config file
rule nextflow_config:
    output:
        "nextflow.config"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 512,
        time_min=lambda wildcards, attempt: attempt * 5,
        tmpdir="tmp"
    params:
        configpath=workflow.source_path("../config/nextflow.config"),
        rsync="--checksum --verbose --human --partial --progress"
    log:
        "logs/nf-core/config.log"
    shell:
        "rsync {params.rsync} {params.configpath} {output} > {log} 2>&1 "