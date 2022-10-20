rule concat_to_bigtable:
    input:
        expand(
            "vep/{annot}/{sample}.tsv",
            annot=["bcr", "hc", "mutect"],
            sample=samples_list,
        ),
    output:
        temp("bigtable/raw.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bigtable/raw.log",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    params:
        "'(NR == 1) || (FNR > 1)'",
    shell:
        "awk {params} {input} > {output} 2> {log}"


rule bigtable_header:
    input:
        "bigtable/raw.tsv",
    output:
        temp("bigtable/header.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
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
        pipe("bigtable/noheader.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bigtable/noheader.log",
    params:
        "'1d'",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sed {params} {input} > {output} 2> {log}"


rule bigtable_sort:
    input:
        "bigtable/noheader.tsv",
    output:
        temp("bigtable/sorted.tsv"),
    threads: 1
    resources:
        mem_mb=get_10gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bigtable/sort.log",
    params:
        "-k1,1 -k2,2n",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "sort {params} {input} > {output}"


rule bigtable_output:
    input:
        header="bigtable/header.tsv",
        content="bigtable/sorted.tsv",
    output:
        protected("data_output/bigtable.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/bigtable/output.log",
    params:
        "",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "cat {input.header} {input.content} > {output} 2> {log}"
