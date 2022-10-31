# Disease Ontology analysis
"""
028.enrichDO
from
-> 026.expand_rank_list
by
-> End job
"""


rule enrichDO:
    input:
        rds="0.26.clusterprofiler/gene_lists/ENTREZID/{comparisons}.RDS",
        universe="026.gene_lists/universe/{comparison}.RDS",
    output:
        rds=temp(
            "027.enrich/DiseaseOnt.ENTREZID/{comparison}/enrich.DiseaseOnt.ENTREZID.RDS"
        ),
        tsv=protected(
            "data_output/GSEA/{comparison}/DiseaseOnt.ENTREZID/enrichment.tsv"
        ),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        enrichDO_extra=config["clusterprofiler"].get(
            "enrich_do", "pvalueCutoff = 1, qvalueCutoff = 1"
        ),
        organism=config.get("organism", "Hs"),
    log:
        "logs/028.enrichdo/{comparison}.log",
    wrapper:
        "bio/clusterProfiler/enrichDO"


# DisGeNET analysis
"""
028.enrichDGN
from
-> 026.expand_rank_list
by
-> End job
"""


rule enrichDGN:
    input:
        rds="026.gene_lists/ENTREZID/{comparison}.RDS",
        universe="026.gene_lists/universe/{comparison}.RDS",
    output:
        rds=temp(
            "027.enrich/DisGenNet.ENTREZID/{comparison}/enrich.DisGenNet.ENTREZID.RDS"
        ),
        tsv=protected("data_output/GSEA/{comparison}/DisGenNet.ENTREZID/enrichment.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        enrichDO_extra=config["clusterprofiler"].get(
            "enrich_dgn", "pvalueCutoff = 1, qvalueCutoff = 1"
        ),
        organism=config.get("organism", "Hs"),
    log:
        "logs/028.enrichdgn/{comparison}.log",
    wrapper:
        "bio/clusterProfiler/enrichDGN"


# Network of Cancer Genes analysis
"""
028.enrichNCG
from
-> 026.expand_rank_list
by
-> End job
"""


rule enrichNCG:
    input:
        rds="026.gene_lists/ENTREZID/{comparison}.RDS",
        universe="026.gene_lists/universe/{comparison}.RDS",
    output:
        rds=temp(
            "027.enrich/NetworkCancerGenes.ENTREZID/{comparison}/enrich.NetworkCancerGenes.ENTREZID.RDS"
        ),
        tsv=protected(
            "data_output/GSEA/{comparison}/NetworkCancerGenes.ENTREZID/enrichment.tsv"
        ),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    params:
        enrichDO_extra=config["clusterprofiler"].get(
            "enrich_ncg", "pvalueCutoff = 1, qvalueCutoff = 1"
        ),
        organism=config.get("organism", "Hs"),
    log:
        "logs/028.enrichncg/{comparison}.log",
    wrapper:
        "bio/clusterProfiler/enrichNCG"
