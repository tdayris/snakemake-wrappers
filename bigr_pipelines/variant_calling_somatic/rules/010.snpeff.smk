snpeff_snpsift_config = {"ref": config["ref"], **config["snpeff_snpsift"]}


module snpeff_meta:
    snakefile:
        str(
            worflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "snpeff_annotate"
            / "test"
            / "Snakefile"
        )
    config:
        snpeff_snpsift_config


use rule snpeff from snpeff_meta with:
    input:
        calls="bcftools/mutect2/{sample}.vcf.gz",
        calls_index=get_tbi("bcftools/mutect2/{sample}.vcf.gz"),
        db=config["ref"]["snpeff"],
