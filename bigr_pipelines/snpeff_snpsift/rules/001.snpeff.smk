# Annotate with snpeff
snpeff_meta_config = {
    "ref": config["ref"],
    **config["snpeff_snpsift"]
}

module snpeff_meta:
    snakefile: "../../meta/bio/snpeff_annotate/test/Snakefile"
    config: snpeff_snpsift_config

use rule snpeff from snpeff_meta with:
    input:
        calls="data_input/calls/{sample}.vcf.gz",
        calls_index="data_input/calls/{sample}.vcf.gz.tbi",
        db=config["ref"]["snpeff"]