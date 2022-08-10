"""
sambamba view 0.6.5 : sambamba view -S -h -f bam -t {threads} -o {sample}.cutadapt.bam {sample}.cutadapt.sam
"""


rule sambamba_view:
    input:
        "bwa/mem/{sample}.sam",
    output:
        temp("sambamba/bam/{sample}.bam"),
    threads: 10
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    conda:
        "envs/sambamba.yaml"
    params:
        "-S -h -f bam",
    log:
        "sambamba/view/{sample}.log",
    shell:
        "sambamba view {params} -t {threads} -o {output} {input} > {log} 2>&1"


"""
sambamba sort 0.6.5 : sambamba sort -t {threads} -o {sample}.cutadapt.sorted.bam {sample}.cutadapt.bam
"""


rule sambamba_sort:
    input:
        "sambamba/bam/{sample}.bam",
    output:
        temp("sambamba/sort/{sample}.bam"),
    threads: 4
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir="tmp",
    conda:
        "envs/sambamba.yaml"
    params:
        "",
    log:
        "logs/sambamba/sort/{sample}.log",
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input} > {log} 2>&1"


"""
sambamba markdup 0.6.5 (optional) : sambamba markdup -r -t {threads} {sample}.cutadapt.sorted.bam {sample}.cutadapt.sorted.rmmarkdup.bam
"""


rule sambamba_markdup:
    input:
        "sambamba/sort/{sample}.bam",
    output:
        "sambamba/markdup/{sample}.bam",
    threads: 10
    resources:
        time_min=lambda wildcards, attempt: attempt * 45,
        mem_mb=lambda wildcards, attempt: attempt * 1024 * 4,
        tmpdir="tmp",
    conda:
        "envs/sambamba.yaml"
    params:
        "-r",
    log:
        "logs/sambamba/markdup/{sample}.log",
    shell:
        "sambamba markdup {params} -t {threads} {input} {output} > {log} 2>&1"
