somatic_tmb_config = {
    "min_coverage": config["tmb"].get("min_coverage"),
    "tmb_highness_threshold": config["tmb"].get("tmb_highness_threshold", 10),
    "allele_depth_keyname": config["tmb"].get("allele_depth_keyname", "AD"),
    "bed": config["ref"]["capture_kit_bed"],
    "sample_list": design["Sample_id"],
}


module somatic_tmb:
    snakefile:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "somatic_tmb"
            / "test"
            / "Snakefile"
        )
    config:
        somatic_tmb_config


use rule * from somatic_tmb


use rule extract_somatic_mutations from somatic_tmb with:
    input:
        vcf="snpeff_snpsift/snpsift/fixed/{sample}.vcf.gz",
