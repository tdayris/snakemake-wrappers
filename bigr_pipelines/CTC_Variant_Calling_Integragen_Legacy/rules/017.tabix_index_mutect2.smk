"""
No command line provided
"""

rule tabix_mutect2_ctc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf.gz",
    output:
        "gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/tabix/{sample}_{version}_{manip}_{nb}.mutect2.log",
    params:
        "-p vcf",
    wrapper:
        "bio/tabix/index"


rule tabix_mutect2_wbc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}.vcf.gz",
    output:
        "gatk/mutect2/{sample}_{version}_{manip}.vcf.gz.tbi",
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/tabix/{sample}_{version}_{manip}.mutect2.log",
    params:
        "-p vcf",
    wrapper:
        "bio/tabix/index"


rule unzip_mutect2_ctc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf.gz",
    output:
        temp("gatk/mutect2/{sample}_{version}_{manip}_{nb}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/gatk/mutect2/unzip/{sample}_{version}_{manip}_{nb}.log",
    params:
        "--decompress --force --verbose --stdout",
    shell:
        "gzip {params} {input} > {output} 2> {log}"


rule unzip_mutect2_wbc:
    input:
        "gatk/mutect2/{sample}_{version}_{manip}.vcf.gz",
    output:
        temp("gatk/mutect2/{sample}_{version}_{manip}.vcf"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir=tmp,
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/gatk/mutect2/unzip/{sample}_{version}_{manip}.log",
    params:
        "--decompress --force --verbose --stdout",
    shell:
        "gzip {params} {input} > {output} 2> {log}"