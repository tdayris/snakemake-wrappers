"""
sambamba markdup 0.6.5 (optional) : sambamba markdup -r -t {threads} {sample}.cutadapt.sorted.bam {sample}.cutadapt.sorted.rmmarkdup.bam
"""


rule sambamba_markdup_ctc:
    input:
        "sambamba/sort/{sample}_{version}_{manip}_{nb}.bam",
    output:
        "sambamba/markdup/{sample}_{version}_{manip}_{nb}.bam",
    threads: 10
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-r",
    log:
        "logs/sambamba/markdup/{sample}_{version}_{manip}_{nb}.log",
    shell:
        "sambamba markdup {params} -t {threads} {input} {output} > {log} 2>&1"


rule sambamba_markdup_wbc:
    input:
        "sambamba/sort/{sample}_{version}_{manip}.bam",
    output:
        "sambamba/markdup/{sample}_{version}_{manip}.bam",
    threads: 10
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-r",
    log:
        "logs/sambamba/markdup/{sample}_{version}_{manip}.log",
    shell:
        "sambamba markdup {params} -t {threads} {input} {output} > {log} 2>&1"


rule sambamba_markdup_baseline:
    input:
        "sambamba/sort/{sample}.baseline.bam",
    output:
        "sambamba/markdup/{sample}.baseline.bam",
    threads: 10
    resources:
        time_min=get_45min_per_attempt,
        mem_mb=get_8gb_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-r",
    log:
        "logs/sambamba/markdup/{sample}.baseline.log",
    shell:
        "sambamba markdup {params} -t {threads} {input} {output} > {log} 2>&1"