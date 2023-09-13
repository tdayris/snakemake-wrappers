"""
cutadapt (v1.14) : cutadapt -a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -o {sample}.cutadapt.fastq -e 0.2 {sample}.fastq
"""


rule cutadapt_ctc:
    input:
        "fastq/{sample}_{version}_{manip}_{nb}.fastq",
    output:
        temp("cutadapt/{sample}_{version}_{manip}_{nb}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    conda:
        str(workflow_source_dir / "envs" / "cutadapt.yaml")
    params:
        "-a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -e 0.2",
    log:
        "logs/cutadapt/{sample}_{version}_{manip}_{nb}.log",
    shell:
        "cutadapt {params} -o {output} {input} > {log} 2>&1"


rule cutadapt_wbc:
    input:
        "fastq/{sample}_{version}_{manip}.fastq",
    output:
        temp("cutadapt/{sample}_{version}_{manip}.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    conda:
        str(workflow_source_dir / "envs" / "cutadapt.yaml")
    params:
        "-a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -e 0.2",
    log:
        "logs/cutadapt/{sample}_{version}_{manip}.log",
    shell:
        "cutadapt {params} -o {output} {input} > {log} 2>&1"



rule cutadapt_baseline:
    input:
        "fastq/{sample}.baseline.fastq",
    output:
        temp("cutadapt/{sample}.baseline.fastq"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    group:
        "clean_input"
    conda:
        str(workflow_source_dir / "envs" / "cutadapt.yaml")
    params:
        "-a ACTGACAGCAGGAATCCCACT -g AGTGGGATTCCTGCTGTCAGT -n 2 -e 0.2",
    log:
        "logs/cutadapt/{sample}.baseline.log",
    shell:
        "cutadapt {params} -o {output} {input} > {log} 2>&1"