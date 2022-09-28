# This rule imports counts from tables to R data object. Its memory requirements
# are linked to the number of samples
"""
007.tximport
from
-> 004.tx_to_gene
-> 004.salmon_quant
by
-> 008.deseq2_dataset_from_tximport
"""


rule tximport:
    input:
        quant=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/quant.sf",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        quant_genes=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/quant.genes.sf",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        lib=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/lib_format_counts.json",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        aux_info=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/aux_info",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        cmd_info=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/cmd_info.json",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        libparams=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/libParams",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        logs=lambda wildcards: expand(
            "004.salmon/pseudo_mapping/{sample}/logs",
            sample=samples_per_prefixes[wildcards.comparison],
        ),
        tx_to_gene="004.salmon/tx2gene_small.tsv",
    output:
        txi=temp("007.tximport/txi.{comparison}.RDS"),
    threads: 1
    resources:
        mem_mb=get_15gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["deseq2"].get(
            "tximport_extra",
            "type='salmon', ignoreTxVersion=TRUE, ignoreAfterBar=TRUE",
        ),
        sample_names=lambda wildcards: samples_per_prefixes[wildcards.comparison],
    log:
        "logs/007.tximport/{comparison}.log",
    wrapper:
        "bio/tximport"
