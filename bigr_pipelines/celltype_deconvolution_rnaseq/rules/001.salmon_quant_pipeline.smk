snakefile_salmon_quant = snakemake.workflow.srcdir(
    "../salmon_quant/Snakefile"
)

module salmon_quant_workflow:
    snakefile: snakefile_salmon_quant
    config: config


use rule * from salmon_quant_workflow