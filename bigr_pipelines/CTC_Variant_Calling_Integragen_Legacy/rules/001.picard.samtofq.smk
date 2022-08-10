"""
picard (v2.23.8) (optional depending on the provided filetype) : picard SamToFastq --INPUT {input}.bam --FASTQ {input}.fastq
"""


rule picard_sam_to_fastq:
    input:
        "data_input/{sample}.{status}.bam",
    output:
        temp("fastq/{sample}.{status}.fastq"),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1024 * attempt,
        time_min=lambda wildcards, attempt: 10 * attempt,
        tmpdir="tmp",
    group:
        "clean_input"
    log:
        "logs/picard/samtofq/{sample}.{status}.log",
    conda:
        "envs/picard.yaml"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output} > {log} 2>&1"
