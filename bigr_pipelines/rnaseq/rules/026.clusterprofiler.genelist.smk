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
        tsv="008.deseq2/{comparison}/wald.{comparison}.tsv",
    output:
        tsv=temp("026.clusterprofiler/{comparison}/wald.{comparison}.tsv"),
    threads: 1
    resources:
        mem_mb=get_768mb_per_attempt,
        time_min=get_10min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/026.clusterprofiler/make_rank_list.{comparison}.log",
    params:
        "-f1,3",
    conda:
        str(workflow_source_dir / "envs" / "bash.yaml")
    shell:
        "cut {params} {input.tsv} > {output.tsv} 2> {log}"


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
        tsv="026.clusterprofiler/{comparison}/wald.{comparison}.tsv",
    output:
        tsv_ensembl=temp("026.clusterprofiler/gene_lists/ENSEMBL/{comparison}.tsv"),
        tsv_ensemblprot=temp(
            "026.clusterprofiler/gene_lists/ENSEMBLPROT/{comparison}.tsv"
        ),
        tsv_symbol=temp("026.clusterprofiler/gene_lists/SYMBOL/{comparison}.tsv"),
        tsv_entrez=temp("026.clusterprofiler/gene_lists/ENTREZID/{comparison}.tsv"),
        rds_ensembl=temp("026.clusterprofiler/gene_lists/ENSEMBL/{comparison}.RDS"),
        rds_ensemblprot=temp(
            "026.clusterprofiler/gene_lists/ENSEMBLPROT/{comparison}.RDS"
        ),
        rds_symbol=temp("026.clusterprofiler/gene_lists/SYMBOL/{comparison}.RDS"),
        rds_entrez=temp("026.clusterprofiler/gene_lists/ENTREZID/{comparison}.RDS"),
        universe_entrez=temp("026.clusterprofiler/universe/ENTREZID/{comparison}.tsv"),
        universe_ensembl=temp("026.clusterprofiler/universe/ENSEMBL/{comparison}.tsv"),
        universe_ensemblprot=temp(
            "026.clusterprofiler/universe/ENSEMBLPROT/{comparison}.tsv"
        ),
        universe_symbol=temp("026.clusterprofiler/universe/SYMBOL/{comparison}.tsv"),
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
        "logs/026.clusterprofiler/expand.{comparison}.log",
    conda:
        str(workflow_source_dir / "envs" / "clusterprofiler.yaml")
    script:
        str(workflow_source_dir / "scripts" / "026.bitr.R")
