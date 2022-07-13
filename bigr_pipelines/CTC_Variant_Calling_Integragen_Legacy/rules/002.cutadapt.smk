"""
cutadapt (v1.14) : cutadapt -a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -o {sample}.cutadapt.fastq -e 0.2 {sample}.fastq
"""
rule cutadapt:
    input:
        "fastq/{sample}.{status}.fastq"
    output:
        temp("cutadapt/{sample}.{status}.fastq")
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 2,
        time_min = lambda wildcards, attempt: attempt * 20,
        tmpdir="tmp"
    group:
        "clean_input"
    conda:
        "envs/cutadapt.yaml"
    params:
        "-a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -e 0.2"
    log:
        "logs/cutadapt/{sample}.{status}.log"
    shell:
        "cutadapt {params} -o {outpuy} {input} > {log} 2>&1"