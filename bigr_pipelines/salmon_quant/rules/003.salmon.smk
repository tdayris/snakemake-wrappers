##############################################
### Configure and call Salmon meta-wrapper ###
##############################################


salmon_config = {
    "genome": config["ref"]["genome"],
    "transcriptome": config["ref"]["transcriptome"],
    "gtf": config["ref"]["gtf"],
    "salmon_libtype": config["params"]["salmon_libtype"],
    "salmon_quant_extra": config["params"]["salmon_quant_extra"],
    "salmon_index_extra": config["params"]["salmon_index_extra"],
    "gentrome": config["ref"]["gentrome"],
    "decoy": config["ref"]["decoys"],
    "index": config["ref"]["index"],
}


module salmon_meta:
    snakefile:
        "../../../meta/bio/salmon/test/Snakefile"
    config:
        salmon_config


use rule * from salmon_meta


use rule salmon_quant_paired from salmon_meta with:
    output:
        quant="salmon/pseudo_mapping/{sample}/quant.sf",
        quant_genes="salmon/pseudo_mapping/{sample}/quant.genes.sf",
        lib="salmon/pseudo_mapping/{sample}/lib_format_counts.json",
