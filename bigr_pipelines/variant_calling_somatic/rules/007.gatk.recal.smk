
gatk_bqsr_config = {
    "threads": config["threads"],
    "genome": config["ref"]["fasta"],
    "dbsnp": config["ref"]["dbsnp"],
    "base_recal_extra": "",
    "apply_base_recal_extra": config.get(
        "gatk", {"apply_base_recal_extra": "--create-output-bam-index"}
    ).get("apply_base_recal_extra", "--create-output-bam-index")
}

module gatk_bqsr_meta:
    snakefile: str(worflow_source_dir / ".." / ".." / "meta" / "bio" / "gatk_bqsr" / "test" / "Snakefile")
    config: gatk_bqsr_config


use rule gatk_apply_baserecalibrator from gatk_bqsr_meta with:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_{status}.bam"),
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
        recal_table="gatk/recal_data_table/{sample}_{status}.grp"
    output:
        bam="gatk/recal_bam/{sample}_{status}.bam",
        bai=get_bai("gatk/recal_bam/{sample}_{status}.bam")
    message:
        "Applying BQSR on {wildcards.status} {wildcards.sample} with GATK"
    params:
        extra=config.get(
            "gatk", {"apply_base_recal_extra", "--create-output-bam-index"}
        ).get("apply_base_recal_extra", "--create-output-bam-index")
    log:
        "logs/gatk/applybqsr/{sample}.{status}.log"


use rule gatk_compute_baserecalibration_table from gatk_bqsr_meta with:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_{status}.bam"),
        ref=config['ref']['fasta'],
        ref_idx=get_fai(config['ref']['fasta']),
        ref_dict=get_dict(config['ref']['fasta']),
        known=config['ref']['dbsnp'],
        known_idx=get_tbi(config['ref']['dbsnp'])
    output:
        recal_table=temp("gatk/recal_data_table/{sample}_{status}.grp")
    message:
        "Compute BQSR table from {wildcards.status} {wildcards.sample} "
        "with GATK"
    log:
        "logs/gatk3/compute_bqsr/{sample}.{status}.log"