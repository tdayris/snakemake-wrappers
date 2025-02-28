"""
sambamba view 0.6.5 : sambamba view -S -h -f bam -t {threads} -o {sample}.cutadapt.bam {sample}.cutadapt.sam
"""


rule sambamba_view_ctc:
    input:
        "bwa/mem/{sample}_{version}_{manip}_{nb}.sam",
    output:
        temp("sambamba/bam/{sample}_{version}_{manip}_{nb}.bam"),
    threads: 10
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-S -h -f bam",
    log:
        "sambamba/view/{sample}_{version}_{manip}_{nb}.log",
    shell:
        "sambamba view {params} -t {threads} -o {output} {input} > {log} 2>&1"


rule sambamba_view_wbc:
    input:
        "bwa/mem/{sample}_{version}_{manip}.sam",
    output:
        temp("sambamba/bam/{sample}_{version}_{manip}.bam"),
    threads: 10
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-S -h -f bam",
    log:
        "sambamba/view/{sample}_{version}_{manip}.log",
    shell:
        "sambamba view {params} -t {threads} -o {output} {input} > {log} 2>&1"


rule sambamba_view_baseline:
    input:
        "bwa/mem/{sample}.baseline.sam",
    output:
        temp("sambamba/bam/{sample}.baseline.bam"),
    threads: 10
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "sambamba.yaml")
    params:
        "-S -h -f bam",
    log:
        "sambamba/view/{sample}.baseline.log",
    shell:
        "sambamba view {params} -t {threads} -o {output} {input} > {log} 2>&1"

