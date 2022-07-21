"""
This rule splits a complex design into small ones.
These simple designs contains only the required factors
and sample names.
"""


rule split_design:
    input:
        design=config["design"],
    output:
        expand("deseq2/designs/{comparison}.tsv", comparison=output_prefixes),
    message:
        "Expanding design in order to make results more readeble"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/deseq2/split_design.log",
    params:
        columns_to_aggregate=config["deseq2"]["design"].get("columns_to_aggregate"),
        columns_to_remove=config["deseq2"]["design"].get("columns_to_remove"),
    wrapper:
        "bio/BiGR/split_design"


"""
This rule formats counts for DESeq2. The design matrix and its corresponding
formula are included.
"""


rule deseq2_dataset_from_tximport:
    input:
        tximport="tximport/txi.{comparison}.RDS",
        coldata="deseq2/designs/{comparison}.tsv",
    output:
        dds=temp("deseq2/{comparison}/dds.{comparison}.RDS"),
    message:
        "Formatting {wildcards.comparison} counts for DESeq2"
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    params:
        design=lambda wildcards: (
            f"~{contrasts[wildcards.comparison][0]}"
            if (batch_effect is False) or wildcards.comparison == "BatchEffect"
            else f"~BatchEffect+{contrasts[wildcards.comparison][0]}"
        ),
        levels=lambda wildcards: contrasts[wildcards.comparison][1:],
        factor=lambda wildcards: contrasts[wildcards.comparison][0],
        ref_level=lambda wildcards: contrasts[wildcards.comparison][-1],
        remove_zeros=True,
        count_filter=5,
    log:
        "logs/deseq2/deseq2_dataset_from_tximport/{comparison}.log",
    wrapper:
        "bio/deseq2/DESeqDataSetFromTximport"


"""
This rule performs the size factor and dispersions estimations as well as the
wald test.
"""


rule deseq2:
    input:
        dds="deseq2/{comparison}/dds.{comparison}.RDS",
    output:
        rds=temp("deseq2/{comparison}/wald.{comparison}.RDS"),
        deseq2_tsv=temp("deseq2/{comparison}/wald.{comparison}.tsv"),
        normalized_counts=temp("deseq2/{comparison}/dst.{comparison}.tsv"),
        dst=temp("deseq2/{comparison}/dst.{comparison}.RDS"),
        intermediar_values=temp("deseq2/{comparison}/mcols.{comparison}.tsv"),
        assays_mu=temp("deseq2/{comparison}/assays.mu.{comparison}.tsv"),
        filter_theta=temp("deseq2/{comparison}/filter.theta.{comparison}.tsv"),
        metadata=temp("deseq2/{comparison}/metadata.{comparison}.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 3
    params:
        contrast=lambda wildcards: contrasts[wildcards.comparison],
    log:
        "logs/deseq2/deseq/{comparison}.log",
    wrapper:
        "bio/deseq2/DESeq"
