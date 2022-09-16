rule term2name_GMT:
    input:
        lambda wildcards: config["gmt"][wildcards.database],
    output:
        temp("clusterprofiler/gmt/{database}.{keytype}.term2name.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/term2name/ensembl/{database}.{keytype}.log"
    params:
        extra=lambda wildcards: get_t2n_extra(wildcards.database)
    shell:
        "awk '{params}' {input} > {output} 2> {log}"


rule enrich_GMT:
    input:
        rds = "gene_lists/{keytype}/{comparison}.RDS",
        universe = "gene_lists/universe/{comparison}.RDS",
        gmt=lambda wildcards: config["gmt"][wildcards.database],
        term2name="clusterprofiler/gmt/{database}.term2name.tsv",
    output:
        readable_rds=temp("enrich/{database}/{comparison}/enrich.{database}.{comparison}.{keytype}.RDS"),
        readable_tsv=protected("data_output/{comparison}/{database}.{keytype}/enrich.{comparison}.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["clusterprofiler"].get("enricher", "pvalueCutoff = 1, qvalueCutoff = 1"),
        org=config.get("organism", "Hs"),
        keytype=lambda wildcards: str(wildcards.keytype)
    log:
        "logs/enricher/{database}/{comparison}.{keytype}.log"
    wrapper:
        "bio/clusterProfiler/enrichGMT"