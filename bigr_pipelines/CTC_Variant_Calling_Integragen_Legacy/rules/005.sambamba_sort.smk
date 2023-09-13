
"""
sambamba sort 0.6.5 : sambamba sort -t {threads} -o {sample}.cutadapt.sorted.bam {sample}.cutadapt.bam
"""


rule sambamba_sort_ctc:
    input:
        "sambamba/bam/{sample}_{version}_{manip}_{nb}.bam",
    output:
        temp("sambamba/sort/{sample}_{version}_{manip}_{nb}.bam"),
    threads: 4
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "",
    log:
        "logs/sambamba/sort/{sample}_{version}_{manip}_{nb}.log",
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input} > {log} 2>&1"



rule sambamba_sort_wbc:
    input:
        "sambamba/bam/{sample}_{version}_{manip}.bam",
    output:
        temp("sambamba/sort/{sample}_{version}_{manip}.bam"),
    threads: 4
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "",
    log:
        "logs/sambamba/sort/{sample}_{version}_{manip}.log",
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input} > {log} 2>&1"


rule sambamba_sort_baseline:
    input:
        "sambamba/bam/{sample}.baseline.bam",
    output:
        temp("sambamba/sort/{sample}.baseline.bam"),
    threads: 4
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "",
    log:
        "logs/sambamba/sort/{sample}.baseline.log",
    shell:
        "sambamba sort {params} -t {threads} -o {output} {input} > {log} 2>&1"