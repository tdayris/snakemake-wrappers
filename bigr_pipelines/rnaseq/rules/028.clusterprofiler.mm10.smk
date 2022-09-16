rule expand_rank_list:
    input:
        tsv=config.get("ranks.list.tsv", "ranks.list.tsv")
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
        time_min=get_10min_per_attempt,
        tmpdir="tmp"
    params:
        gene_id_type=config.get(
            "gene_id_type", "ENSEMBL"
        )
    log:
        "logs/expand.log"
    wrapper:
        "bio/clusterProfiler/mm10_genelist"