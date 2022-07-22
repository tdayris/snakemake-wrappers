vcf_post_process_config = {
    "ncbi_build": config["params"].get("ncbi_build", "GRCh38"),
    "center": config["params"].get("center", "GustaveRoussy"),
    "annotation_tag": "ANN=",
    "sample_list": design["Sample_id"].to_list(),
    "genome": config["ref"]["fasta"],
    "known": config["ref"]["dbsnp"],
    "chr": config["params"]["chr"],
}


module vcf_post_process:
    snakefile:
        str(
            workflow_source_dir
            / ".."
            / ".."
            / "meta"
            / "bio"
            / "vcf_post_process"
            / "test"
            / "Snakefile"
        )
    config:
        vcf_post_process_config


use rule * from vcf_post_process
