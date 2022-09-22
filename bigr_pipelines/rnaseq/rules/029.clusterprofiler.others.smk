"""
Disease Ontology analysis
"""
rule enrichDO:
    input:
        rds = "gene_lists/ENTREZID/{comparison}.RDS",
        universe = "gene_lists/universe/{comparison}.RDS"
    output:
        rds = temp(
            "enrich/DiseaseOnt/{comparison}/enrich.DiseaseOnt.{comparison}.ENTREZID.RDS"
        ),
        tsv = "results/{comparison}/DiseaseOnt.ENTREZID/enrich.{comparison}.tsv"
    message: "Running enrichDO on {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    params:
        enrichDO_extra=config["clusterprofiler"].get("enrich_do", "pvalueCutoff = 1, qvalueCutoff = 1"),
        organism = config.get("organism", "Hs"),
    log:
        "logs/enrichdo/{comparison}.log"
    wrapper:
        "bio/clusterProfiler/enrichDO"


"""
DisGeNET analysis
"""
rule enrichDGN:
    input:
        rds = "gene_lists/ENTREZID/{comparison}.RDS",
        universe = "gene_lists/universe/{comparison}.RDS"
    output:
        rds = temp(
            "enrich/DisGenNet/{comparison}/enrich.DisGenNet.{comparison}.ENTREZID.RDS"
        ),
        tsv = "results/{comparison}/DisGenNet.ENTREZID/enrich.{comparison}.tsv"
    message: "Running enrichDGN on {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    params:
        enrichDO_extra=config["clusterprofiler"].get("enrich_dgn", "pvalueCutoff = 1, qvalueCutoff = 1"),
        organism = config.get("organism", "Hs"),
    log:
        "logs/enrichdgn/{comparison}.log"
    wrapper:
        "bio/clusterProfiler/enrichDGN"


"""
Network of Cancer Genes analysis
"""
rule enrichNCG:
    input:
        rds = "gene_lists/ENTREZID/{comparison}.RDS",
        universe = "gene_lists/universe/{comparison}.RDS"
    output:
        rds = temp(
            "enrich/NetworkCancerGenes/{comparison}/enrich.NetworkCancerGenes.{comparison}.ENTREZID.RDS"
        ),
        tsv = "results/{comparison}/NetworkCancerGenes.ENTREZID/enrich.{comparison}.tsv"
    message: "Running enrichNCG on {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    params:
        enrichDO_extra=config["clusterprofiler"].get("enrich_ncg", "pvalueCutoff = 1, qvalueCutoff = 1"),
        organism = config.get("organism", "Hs"),
    log:
        "logs/enrichncg/{comparison}.log"
    wrapper:
        "bio/clusterProfiler/enrichNCG"