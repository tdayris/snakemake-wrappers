# Expand terms in human readable format for GMT files
"""
027.term2name_GMT
from
-> Entry Job
by
-> 027.enrich_GMT
"""


rule term2name_GMT:
    input:
        lambda wildcards: config["clusterprofiler"]["gmt"][wildcards.database],
    output:
        temp("026.clusterprofiler/gmt/{database}.{keytype}.term2name.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    log:
        "logs/027.term2name/ensembl/{database}.{keytype}.log",
    params:
        extra=lambda wildcards: get_t2n_extra(wildcards.database),
    shell:
        "awk '{params}' {input} > {output} 2> {log}"


# Perform gene set enrichment on GMT files
"""
027.term2name_GMT
from
-> 027.term2name_GMT
-> 026.expand_rank_list
by
-> End job
"""


rule enricher_GMT:
    input:
        gene="026.clusterprofiler/gene_lists/{keytype}/{comparison}.tsv",
        term2gene=lambda wildcards: config["clusterprofiler"]["gmt"][wildcards.database],
        term2name="026.clusterprofiler/gmt/{database}.{keytype}.term2name.tsv",
    output:
        rds=temp(
            "027.enrich/{database}.{keytype}/{comparison}/enrich.{database}.{keytype}.RDS"
        ),
        tsv=protected(
            "data_output/GSEA/{comparison}/{database}.{keytype}/enrichment.tsv"
        ),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["clusterprofiler"].get(
            "enricher", "pvalueCutoff = 1, qvalueCutoff = 1"
        ),
    log:
        "logs/027.enricher/{database}/{comparison}.{keytype}.log",
    wrapper:
        "bio/clusterprofiler/enricher"
