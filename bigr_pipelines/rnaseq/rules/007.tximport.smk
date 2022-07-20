"""
This rule imports counts from tables to R data object. Its memory requirements
are linked to the number of samples
"""
rule tximport:
    input:
        quant=lambda wildcards: expand(
            "salmon/pseudo_mapping/{sample}/quant.sf",
            sample=samples_per_prefixes[wildcards.comparison]
        ),
        tx_to_gene="salmon/tx2gene.tsv"
    output:
        txi=temp("tximport/txi.{comparison}.RDS")
    message: "Importing counts in DESeq2 for {wildcards.comparison}"
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp"
    params:
        extra=config["deseq2"].get(
            "tximport_extra",
            "type='salmon', ignoreTxVersion=TRUE, ignoreAfterBar=TRUE"
        )
    log:
        "logs/tximport/{comparison}.log"
    wrapper:
        "bio/tximport"