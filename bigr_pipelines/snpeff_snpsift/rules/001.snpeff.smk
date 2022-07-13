# Annotate with snpeff
snpeff_meta_config = {
    "ref": config["ref"],
    **config["snpeff_snpsift"]
}

module snpsift:
    snakefile: "../../meta/bio/snpsift/test/Snakefile"
    config: snpeff_meta_config


use rule * from snpsift