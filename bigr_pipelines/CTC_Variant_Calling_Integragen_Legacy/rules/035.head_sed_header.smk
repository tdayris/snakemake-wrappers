"""
No command line provided
"""

rule bigtable_header:
    input:
        "bigtable/raw.tsv",
    output:
        temp("bigtable/header.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/header.log",
    params:
        "-n 1",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "head {params} {input} > {output} 2> {log}"


rule bigtable_noheader:
    input:
        "bigtable/raw.tsv",
    output:
        temp("bigtable/noheader.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/noheader.log",
    params:
        "'1d;s|.ctc.brc.vcf||g;s|.ctc.bcr.vcf||g;s|Baseline|Germline|g'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params} {input} > {output} 2> {log}"