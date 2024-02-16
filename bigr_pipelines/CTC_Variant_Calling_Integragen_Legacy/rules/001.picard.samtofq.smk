"""
picard (v2.23.8) (optional depending on the provided filetype) : picard SamToFastq --INPUT {input}.bam --FASTQ {input}.fastq
"""


rule picard_sam_to_fastq_ctc:
    input:
        lambda wildcards: get_ctc(wildcards),
    output:
        temp("fastq/{sample}_{version}_{manip}_{nb}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    log:
        "logs/picard/samtofq/{sample}_{version}_{manip}_{nb}.log",
    conda:
        str(workflow_source_dir / "envs" / "picard.yaml")
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output} > {log} 2>&1"


rule picard_sam_to_fastq_wbc:
    input:
        lambda wildcards: get_wbc(wildcards),
    output:
        temp("fastq/{sample}_{version}_{manip}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    log:
        "logs/picard/samtofq/{sample}_{version}_{manip}.log",
    conda:
        str(workflow_source_dir / "envs" / "picard.yaml")
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output} > {log} 2>&1"


rule picard_sam_to_fastq_baseline:
    input:
        lambda wildcards: get_baseline(wildcards),
    output:
        temp("fastq/{sample}.baseline.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    log:
        "logs/picard/samtofq/{sample}.baseline.log",
    conda:
        str(workflow_source_dir / "envs" / "picard.yaml")
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output} > {log} 2>&1"
