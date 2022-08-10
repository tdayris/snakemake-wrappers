"""
cutadapt (v1.14) : cutadapt -a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -o {sample}.cutadapt.fastq -e 0.2 {sample}.fastq
"""


rule cutadapt:
    input:
        "fastq/{sample}.fastq",
    output:
        temp("cutadapt/{sample}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    group:
        "clean_input"
    conda:
        "envs/cutadapt.yaml"
    params:
        "-a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -e 0.2",
    log:
        "logs/cutadapt/{sample}.{status}.log",
    shell:
        "cutadapt {params} -o {output} {input} > {log} 2>&1"
