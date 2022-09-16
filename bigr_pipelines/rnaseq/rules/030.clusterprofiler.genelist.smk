rule make_rank_list:
    input:
        tsv=lambda wildcards: expand(
            "deseq2/{comparison}/wald.{comparison}.tsv",
            comparison=output_prefixes
        ),
    output:
        temp("clusterprofiler/rank.list.tsv")
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/clusterprofiler/make_rank_list.log"
    params:
        rank_on=config["clusterprofiler"].get("rank_on", "padj")
    conda:
        str(worflow_source_dir / "envs" / "bash.yaml")
    script:
        str(worflow_source_dir / "scripts" / "030.make_rank_list.py")


rule expand_rank_list:
    input:
        tsv="clusterprofiler/rank.list.tsv"
    output:
        tsv=expand(
            "gene_lists/ENTREZID/{comparison}.tsv",
            comparison=comparisons
        ),
        entrez_rds=temp(expand(
            "gene_lists/ENTREZID/{comparison}.RDS",
            comparison=comparisons
        )),
        symbol_rds=temp(expand(
            "gene_lists/SYMBOL/{comparison}.RDS",
            comparison=comparisons
        )),
        ensembl_rds=temp(expand(
            "gene_lists/ENSEMBL/{comparison}.RDS",
            comparison=comparisons
        )),
        protein_rds = temp(expand(
            "gene_lists/ENSEMBLPROT/{comparison}.RDS",
            comparison=comparisons
        )),
        universe = temp(expand(
            "gene_lists/universe/{comparison}.RDS",
            comparison=comparisons
        ))
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp"
    params:
        gene_id_type=config.get(
            "gene_id_type", "ENSEMBL"
        )
    log:
        "logs/clusterprofiler/expand.log"
    wrapper:
        "bio/clusterProfiler/hg38_genelist"