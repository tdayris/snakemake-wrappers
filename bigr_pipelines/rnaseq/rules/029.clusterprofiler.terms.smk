# Format term / gene relation from TSV
"""
029.term2gene_TERMS
from
-> Entry Job
by
-> 029.enricher_TERMS
"""


rule term2gene_TERMS:
    input:
        lambda wildcards: config["ppi"][wildcards.database],
    output:
        temp("026.clusterprofiler/term2gene/{database}.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    params:
        begin='FS=OFS="\t"',
        body=["print $3 FS $1"],
    group:
        "prepare_terms"
    log:
        "logs/029.awk/prepare_terms2gene/{database}.log",
    wrapper:
        "bio/awk"


# Expand terms in human readable format for TSV files
"""
029.terms2name_TERMS
from
-> Entry Job
by
-> 029.enricher_TERMS
"""


rule terms2name_TERMS:
    input:
        lambda wildcards: config["ppi"][wildcards.database],
    output:
        temp("026.clusterprofiler/term2name/{database}.tsv"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    params:
        begin='FS=OFS="\t"',
        body=["print $3 FS $4"],
    group:
        "prepare_terms"
    log:
        "logs/029.awk/prepare_terms2name/{database}.log",
    wrapper:
        "bio/awk"


# Perform term enrichment
"""
029.enricher_TERMS
from
-> 029.terms2name_TERMS
-> 029.term2gene_TERMS
by
-> End job
"""


rule enricher_TERMS:
    input:
        rds="026.gclusterprofiler/ene_lists/ENSEMBLPROT/{comparison}.RDS",
        universe="026.clusterprofiler/gene_lists/universe/{comparison}.RDS",
        term_name="026.clusterprofiler/term2name/{database}.tsv",
        term_gene="026.clusterprofiler/term2gene/{database}.tsv",
    output:
        readable_rds=temp(
            "027.enrich/{database}.ENSEMBLPROT/{comparison}/enrich.{database}.ENSEMBLPROT.RDS"
        ),
        readable_tsv=protected("data_output/GSEA/{comparison}/{database}.ENSEMBLPROT/enrichment.tsv"),
    threads: 1
    resources:
        mem_mb=get_3gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["clusterprofiler"].get(
            "enrich_gmt", "pvalueCutoff = 1, qvalueCutoff = 1"
        ),
        org=config.get("organism", "Hs"),
        keytype="ENSEMBLPROT",
    log:
        "logs/029.enricher/{database}/{comparison}.ENSEMBLPROT.log",
    wrapper:
        "bio/clusterProfiler/enrichGMT"
