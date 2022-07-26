"""
This rule indexes the recalibrated bam
"""


rule samtools_index_gatk:
    input:
        "gatk/recal_bam/{sample}.bam",
    output:
        temp("gatk/recal_bam/{sample}.bam.bai"),
    threads: min(config.get("max_threads", 20), 4)
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    params:
        extra=config["samtools"].get("index_extra", ""),
    log:
        "logs/samtools/index/{sample}.gatk.log",
    wrapper:
        "bio/samtools/index"


"""
This rule applies the BQSR to the mapped reads
"""


rule gatk_apply_baserecalibrator:
    input:
        bam="star/{sample}/variants/{sample}.bam",
        bam_index="star/{sample}/variants/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        recal_table="gatk/recal_data_table/{sample}.grp",
    output:
        bam=temp("gatk/recal_bam/{sample}.bam"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk/applybqsr/{sample}.log",
    params:
        extra=config.get("applybqsr", ""),
    wrapper:
        "bio/gatk/applybqsr"


"""
This rule computes BQSR on mapped reads, given a knoledge database
"""


rule gatk_compute_baserecalibration_table:
    input:
        bam="star/{sample}/variants/{sample}.bam",
        bam_index="star/{sample}/variants/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        known=config["reference"]["dbsnp"],
        known_idx=config["refernce"]["dbsnp_tbi"],
    output:
        recal_table=temp("gatk/recal_data_table/{sample}.grp"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/gatk3/compute_bqsr/{sample}.log",
    params:
        extra=config.get("baserecalibrator", ""),
    wrapper:
        "bio/gatk/baserecalibrator"
