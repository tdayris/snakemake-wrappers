"""
sambamba view 0.6.5 : sambamba view -S -h -f bam -t {threads} -o {sample}.cutadapt.bam {sample}.cutadapt.sam
"""
rule sambamba_view:
    input:
        "bwa/mem/{sample}.{status}.sam"
    output:
        temp("sambamba/bam/{sample}.{status}.bam")
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 5,
        time_min = lambda wildcards, attempt: attempt * 25,
        tmpdir = "tmp"
    group:
        "map_n_sort"
    conda:
        "envs/sambamba.yaml"
    params:
        "-S -h -f bam"
    log:
        "sambamba/view/{sample}.{status}.log"
    shell:
        "sambamba view {params} -t {threads} -o {output} {input} > {log} 2>&1"
        
"""
sambamba sort 0.6.5 : sambamba sort -t {threads} -o {sample}.cutadapt.sorted.bam {sample}.cutadapt.bam
"""
rule sambamba_sort:
    input:
        "sambamba/bam/{sample}.{status}.bam"
    output:
        temp("sambamba/sort/{sample}.{status}.bam")
    threads: 4
    resources:
        time_min = lambda wildcards, attempt: attempt * 45,
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 9,
        tmpdir = "tmp"
    group:
        "map_n_sort"
    conda:
        "envs/sambamba.yaml"
    params:
        ""
    log:
        "logs/sambamba/sort/{sample}.{status}.log"
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input} > {log} 2>&1"
        
"""
sambamba markdup 0.6.5 (optional) : sambamba markdup -r -t {threads} {sample}.cutadapt.sorted.bam {sample}.cutadapt.sorted.rmmarkdup.bam
"""
rule sambamba_markdup:
    input:
        "sambamba/sort/{sample}.{status}.bam"
    output:
        "sambamba/markdup/{sample}.{status}.bam"
    threads: 10
    resources:
        time_min = lambda wildcards, attempt: attempt * 45,
        mem_mb = lambda wildcards, attempt: attempt * 1024 * 4
        tmpdir = "tmp"
    conda:
        "envs/sambamba.yaml"
    params:
        "-r"
    log:
        "logs/sambamba/markdup/{sample}.{status}.log"
    shell:
        "sambamba markdup {params} -t {threads} {input} {output} > {log} 2>&1"
    