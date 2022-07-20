"""
This rule splits a complex design into small ones.
These simple designs contains only the required factors
and sample names.
"""
rule split_design:
    input:
        design=config["design"]["path"],
    output:
        expand(
            "deseq2/designs/{comparison}.tsv",
            comparison=output_prefixes
        )
    message:
        "Expanding design in order to make results more readeble"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/deseq2/split_design.log"
    params:
        columns_to_aggregate=config["deseq2"]["design"].get("columns_to_aggregate"),
        columns_to_remove=config["deseq2"]["design"].get("columns_to_remove")
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
        dds=temp("deseq2/{comparison}/dds.{comparison}.RDS")
    message: "Formatting {wildcards.comparison} counts for DESeq2",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 3072, 20480),
        time_min=lambda wildcards, attempt: attempt * 45
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
        count_filter=5
    log:
        "logs/deseq2/deseq2_dataset_from_tximport/{comparison}.log"
    wrapper:
        "bio/deseq2/DESeqDataSetFromTximport"


"""
This rule performs the size factor and dispersions estimations as well as the
wald test.
"""
rule deseq2:
    input:
        dds="deseq2/{comparison}/dds.{comparison}.RDS"
    output:
        rds=temp("deseq2/{comparison}/wald.{comparison}.RDS"),
        deseq2_tsv=report(
            "deseq2/{comparison}/wald.{comparison}.tsv",
            caption=deseq2_rst,
            category="DESeq2 results"
        ),
        normalized_counts=report(
            "deseq2/{comparison}/dst.{comparison}.tsv",
            caption=normalized_counts_rst ,
            category="Normalized counts"
        ),
        dst=temp("deseq2/{comparison}/dst.{comparison}.RDS"),
        intermediar_values="deseq2/{comparison}/mcols.{comparison}.tsv",
        assays_mu="deseq2/{comparison}/assays.mu.{comparison}.tsv",
        filter_theta="deseq2/{comparison}/filter.theta.{comparison}.tsv",
        metadata="deseq2/{comparison}/metadata.{comparison}.tsv"
    message: "Running DESeq2 analysis for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 2048,
        time_min=lambda wildcards, attempt: attempt * 60
    params:
        contrast=lambda wildcards: contrasts[wildcards.comparison],
    log:
        "logs/deseq2/deseq/{comparison}.log"
    wrapper:
        "bio/deseq2/DESeq"