# Build gene rank list from DESeq2
"""
026.make_rank_list
from
-> 008.deseq2
by
-> 026.expand_rank_list
"""


rule make_rank_list:
    input:
        tsv=lambda wildcards: expand(
            "008.deseq2/{comparison}/wald.{comparison}.tsv", comparison=output_prefixes
        ),
    output:
        tsv=temp("026.clusterprofiler/rank.list.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/026.clusterprofiler/make_rank_list.log",
    params:
        rank_on=config["clusterprofiler"].get("rank_on", "padj"),
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    script:
        str(workflow_source_dir / "scripts" / "030.make_rank_list.py")


# Expand and annotate gene rank list
"""
026.expand_rank_list
from
-> 026.make_rank_list
by
-> 027.term2name_GMT
-> 028.enrichDO
-> 028.enrichDGN
-> 028.enrichNCG
-> enricher_TERMS
"""


rule expand_rank_list:
    input:
        tsv="026.clusterprofiler/rank.list.tsv",
    output:
        tsv=expand(
            "026.clusterprofiler/gene_lists/ENTREZID/{comparison}.tsv",
            comparison=output_prefixes,
        ),
        entrez_rds=temp(
            expand(
                "026.clusterprofiler/gene_lists/ENTREZID/{comparison}.RDS",
                comparison=output_prefixes,
            )
        ),
        symbol_rds=temp(
            expand(
                "026.clusterprofiler/gene_lists/SYMBOL/{comparison}.RDS",
                comparison=output_prefixes,
            )
        ),
        ensembl_rds=temp(
            expand(
                "026.clusterprofiler/gene_lists/ENSEMBL/{comparison}.RDS",
                comparison=output_prefixes,
            )
        ),
        protein_rds=temp(
            expand(
                "026.clusterprofiler/gene_lists/ENSEMBLPROT/{comparison}.RDS",
                comparison=output_prefixes,
            )
        ),
        universe=temp(
            expand(
                "026.clusterprofiler/gene_lists/universe/{comparison}.RDS",
                comparison=output_prefixes,
            )
        ),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    params:
        gene_id_type=config.get("gene_id_type", "ENSEMBL"),
        orgdb=(
            "org.Hs.eg.db"
            if config["clusterprofiler"].get("organism", "Hs") == "Hs"
            else "org.Mm.eg.db"
        ),
    log:
        "logs/026.clusterprofiler/expand.log",
    wrapper:
        "bio/clusterProfiler/hg38_genelist"