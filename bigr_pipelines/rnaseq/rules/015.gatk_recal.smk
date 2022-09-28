# This rule indexes the recalibrated bam
"""
015.samtools_index_gatk
from:
-> 015.gatk_apply_baserecalibrator
by:
-> 016.mutect2_germline
"""
rule 015_samtools_index_gatk:
    input:
        "010.gatk/recal_bam/{sample}.bam",
    output:
        temp("010.gatk/recal_bam/{sample}.bam.bai"),
    threads: min(config.get("max_threads", 20), 4)
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["samtools"].get("index_extra", ""),
    log:
        "logs/015.samtools/index/{sample}.gatk.log",
    wrapper:
        "bio/samtools/index"



# This rule applies the BQSR to the mapped reads
"""
015.gatk_apply_baserecalibrator
from:
-> 015.gatk_compute_baserecalibration_table
-> 010.gatk_split_n_cigar_reads
by:
-> 016.mutect2_germline
"""
rule 015_gatk_apply_baserecalibrator:
    input:
        bam="010.gatk/splitncigarreads/{sample}.bam",
        bam_index="010.gatk/splitncigarreads/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        recal_table="010.gatk/recal_data_table/{sample}.grp",
    output:
        bam=temp("010.gatk/recal_bam/{sample}.bam"),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/015.gatk/applybqsr/{sample}.log",
    params:
        extra=config.get("applybqsr", ""),
    wrapper:
        "bio/gatk/applybqsr"



# This rule computes BQSR on mapped reads, given a knoledge database
"""
015.gatk_compute_baserecalibration_table
from:
-> 010.gatk_split_n_cigar_reads
by:
-> 015.gatk_apply_baserecalibrator
"""
rule 015_gatk_compute_baserecalibration_table:
    input:
        bam="010.gatk/splitncigarreads/{sample}.bam",
        bam_index="010.gatk/splitncigarreads/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        known=config["reference"]["dbsnp"],
        known_idx=config["refernce"]["dbsnp_tbi"],
    output:
        recal_table=temp("010.gatk/recal_data_table/{sample}.grp"),
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/015.gatk/compute_bqsr/{sample}.log",
    params:
        extra=config.get("baserecalibrator", ""),
    wrapper:
        "bio/gatk/baserecalibrator"
