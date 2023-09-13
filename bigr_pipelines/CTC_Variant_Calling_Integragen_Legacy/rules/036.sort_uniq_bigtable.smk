"""
No command line provided
"""

rule bigtable_sort:
    input:
        "bigtable/noheader.tsv",
    output:
        temp("bigtable/sorted.tsv"),
    threads: 2
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/sort.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sort {params} {input} | uniq > {output}"



rule bigtable_output:
    input:
        header="bigtable/header.tsv",
        content="bigtable/sorted.tsv",
    output:
        protected("bigtable/duplicated.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/bigtable/output.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "cat {input.header} {input.content} > {output} 2> {log}"



rule remove_duplicates:
    input:
        "data_output/bigtable.tsv"
    output:
        temp("bigtable/uniq.csv")
    threads: 3
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir=tmp,
    log:
        "logs/xsv/preset/bigtable.log"
    params:
        s="",
        u="",
        e="'s|,|;|g;s|\t|,|g'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params.e} {input} | xsv sort {params.s} | uniq {params.u} > {output} 2> {log}"