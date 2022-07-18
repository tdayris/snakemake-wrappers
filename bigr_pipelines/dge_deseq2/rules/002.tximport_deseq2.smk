deseq2_config = {
    "gtf": config["ref"]["gtf"],
    "design": config["design"],
    "output_prefixes": output_prefixes,
    "comparison_levels": comparison_levels,
    "samples_per_prefixes": samples_per_prefixes
}

snakefile_tximport_deseq2 = snakemake.workflow.srcdir(
    "../../meta/bio/tximport_deseq2/test/Snakefile"
)


module tximport_deseq2:
    snakefile: snakefile_tximport_deseq2
    config: deseq2_config


use rule * from tximport_deseq2