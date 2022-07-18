module salmon_quant_workflow:
    snakefile: snakemake.workflow.srcdir("../salmon_quant/Snakefile")
    config: config


use rule * from salmon_quant_workflow as salmon_quant_workflow_*