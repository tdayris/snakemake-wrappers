"""
picard (v2.23.8) (optional depending on the provided filetype) : picard SamToFastq --INPUT {input}.bam --FASTQ {input}.fastq
"""


rule picard_sam_to_fastq:
    input:
        "data_input/{sample}.bam",
    output:
        temp("fastq/{sample}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    group:
        "clean_input"
    log:
        "logs/picard/samtofq/{sample}.log",
    conda:
        "envs/picard.yaml"
    shell:
        "picard SamToFastq --INPUT {input} --FASTQ {output} > {log} 2>&1"
