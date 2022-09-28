# This rule splits a complex design into small ones.
# These simple designs contains only the required factors
# and sample names.
"""
008.split_design
from:
-> Entry Job
by:
-> 008.deseq2_dataset_from_tximport
"""
rule 008_split_design:
    input:
        design=config["design"],
    output:
        expand("008.deseq2/designs/{comparison}.tsv", comparison=output_prefixes),
    message:
        "Expanding design in order to make results more readeble"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/008.deseq2/split_design.log",
    params:
        columns_to_aggregate=config["deseq2"]["design"].get("columns_to_aggregate"),
        columns_to_remove=config["deseq2"]["design"].get("columns_to_remove"),
    wrapper:
        "bio/BiGR/split_design"



# This rule formats counts for DESeq2. The design matrix and its corresponding
# formula are included.
"""
008.deseq2_dataset_from_tximport
from:
-> 008.split_design
"""
rule 008_deseq2_dataset_from_tximport:
    input:
        tximport="007.tximport/txi.{comparison}.RDS",
        coldata="008.deseq2/designs/{comparison}.tsv",
    output:
        dds=temp("008.deseq2/{comparison}/dds.{comparison}.RDS"),
    message:
        "Formatting {wildcards.comparison} counts for DESeq2"
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
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
        count_filter=0.01,
    log:
        "logs/008.deseq2/deseq2_dataset_from_tximport/{comparison}.log",
    wrapper:
        "bio/deseq2/DESeqDataSetFromTximport"


# This rule performs the size factor and dispersions estimations as well as the
# wald test.
"""
008.deseq2
from:
-> 008.deseq2_dataset_from_tximport
by:
-> 009.deseq2_readable
-> 026.make_rank_list
"""
rule 008_deseq2:
    input:
        dds="008.deseq2/{comparison}/dds.{comparison}.RDS",
    output:
        rds=temp("008.deseq2/{comparison}/wald.{comparison}.RDS"),
        deseq2_tsv=temp("008.deseq2/{comparison}/wald.{comparison}.tsv"),
        normalized_counts=temp("008.deseq2/{comparison}/dst.{comparison}.tsv"),
        dst=temp("008.deseq2/{comparison}/dst.{comparison}.RDS"),
        intermediar_values=temp("008.deseq2/{comparison}/mcols.{comparison}.tsv"),
        assays_mu=temp("008.deseq2/{comparison}/assays.mu.{comparison}.tsv"),
        filter_theta=temp("008.deseq2/{comparison}/filter.theta.{comparison}.tsv"),
        metadata=temp("008.deseq2/{comparison}/metadata.{comparison}.tsv"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        contrast=lambda wildcards: contrasts[wildcards.comparison],
    log:
        "logs/008.deseq2/deseq/{comparison}.log",
    wrapper:
        "bio/deseq2/DESeq"
