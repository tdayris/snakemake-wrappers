rule term2name_TERMS:
    input:
        lambda wildcards: config["ppi"][wildcards.database]
    output:
        temp("clusterprofiler/term2gene/{database}.tsv")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp"
    params:
        begin='FS=OFS="\t"',
        body=['print $3 FS $1']
    group:
        "prepare_terms"
    log:
        "logs/awk/prepare_terms2gene/{database}.log"
    wrapper:
        "bio/awk"


rule terms2name_TERMS:
    input:
        lambda wildcards: config["ppi"][wildcards.database]
    output:
        temp("clusterprofiler/term2name/{database}.tsv")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp"
    params:
        begin='FS=OFS="\t"',
        body=['print $3 FS $4']
    group:
        "prepare_terms"
    log:
        "logs/awk/prepare_terms2name/{database}.log"
    wrapper:
        "bio/awk"


rule enricher_TERMS:
    input:
        rds = "gene_lists/ENSEMBLPROT/{comparison}.RDS",
        universe = "gene_lists/universe/{comparison}.RDS",
        term_name = "clusterprofiler/term2name/{database}.tsv",
        term_gene = "clusterprofiler/term2gene/{database}.tsv"
    output:
        readable_rds = temp("enrich/{database}/{comparison}/enrich.{database}.{comparison}.ENSEMBLPROT.RDS"),
        readable_tsv = "results/{comparison}/{database}.ENSEMBLPROT/enrich.{comparison}.ENSEMBLPROT.tsv"
    threads: 1
    resources:
        mem_mb=get_3gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    params:
        extra = config.get("enricher_extra", {"enricher": "pvalueCutoff = 1, qvalueCutoff = 1"}).get("enricher", "pvalueCutoff = 1, qvalueCutoff = 1"),
        org = config.get("organism", "Hs"),
        keytype = "ENSEMBLPROT"
    log:
        "logs/enricher/{database}/{comparison}.ENSEMBLPROT.log"
    wrapper:
        "bio/clusterProfiler/enrichGMT"