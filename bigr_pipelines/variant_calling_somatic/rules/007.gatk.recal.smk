rule gatk_apply_baserecalibrator:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_{status}.bam"),
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        recal_table="gatk/recal_data_table/{sample}_{status}.grp",
    output:
        bam=protected("data_output/BAM/{sample}_{status}.bam"),
        bai=protected("data_output/BAM/{sample}_{status}.bam.bai"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["gatk"].get("apply_base_recal", "--create-output-bam-index"),
    log:
        "logs/gatk/applybqsr/{sample}.{status}.log",
    wrapper:
        "bio/gatk/applybqsr"


rule gatk_compute_baserecalibration_table:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bam_index=get_bai("sambamba/markdup/{sample}_{status}.bam"),
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
        known=config["reference"]["dbsnp"],
        known_idx=config["reference"]["dbsnp_tbi"],
    output:
        recal_table=temp("gatk/recal_data_table/{sample}_{status}.grp"),
    threads: 1
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk3/compute_bqsr/{sample}.{status}.log",
    params:
        extra=config["gatk"].get("base_recalibrator"),
    wrapper:
        "bio/gatk/baserecalibrator"
