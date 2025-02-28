# Acquire multiple stats over BAM files
"""
011.samtools_stats
from
-> 010.sambamba_sort_star
by
-> snakefile.star_fusion_results
"""


rule samtools_stats:
    input:
        aln="010.star/{sample}/{maptype}/{sample}.bam",
        aln_idx="010.star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("011.samtools/stats/{sample}.{maptype}.stats"),
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/011.samtools/stats/{sample}.{maptype}.log",
    params:
        extra=config["samtools"].get("stats", ""),
    wrapper:
        "bio/samtools/stats"


# Index star bam files
"""
011.samtools_index_bam
from
-> 010.sambamba_sort_star
by
-> 011.samtools_stats
"""


rule samtools_index_bam:
    input:
        "star/{sample}/{maptype}/{sample}.bam",
    output:
        temp("star/{sample}/{maptype}/{sample}.bam.bai"),
    threads: min(config.get("max_threads", 20), 4)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["samtools"].get("index_extra", ""),
    log:
        "logs/011.samtools/index/{sample}.{maptype}.log",
    wrapper:
        "bio/samtools/index"


# CRAM star bam files for space saving
"""
011.samtools_cram
from
-> 010.sambamba_sort_star
by
-> End Job
"""


rule samtools_cram:
    input:
        aln="star/{sample}/{maptype}/{sample}.bam",
        aln_idx="star/{sample}/{maptype}/{sample}.bam.bai",
        fasta=config["reference"]["genome"],
        fasta_index=config["reference"]["genome_index"],
    output:
        protected("data_output/Mapping/{maptype}/{sample}.cram"),
    threads: min(config.get("max_threads", 2), 2)
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_4h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["samtools"].get("cram_extra", ""),
    log:
        "logs/011.samtools/cram/{sample}.{maptype}.log",
    wrapper:
        "bio/samtools/view"
