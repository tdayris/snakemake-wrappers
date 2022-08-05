rule samtools_stats:
    input:
        bam="star/{sample}/{maptype}/{sample}.bam",
        bai="star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        bed=config["reference"]["capture_kit_bed"],
    output:
        temp("samtools/stats/{sample}.{maptype}.stats"),
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_35min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/samtools/stats/{sample}.{maptype}.log",
    params:
        extra=config["samtools"].get("stats", ""),
    wrapper:
        "bio/samtools/stats"


rule samtools_index_bam:
    input:
        "star/{sample}/{maptype}/{sample}.bam",
    output:
        temp("star/{sample}/{maptype}/{sample}.bam.bai"),
    threads: min(config.get("max_threads", 20), 4)
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["samtools"].get("index_extra", ""),
    log:
        "logs/samtools/index/{sample}.{maptype}.log",
    wrapper:
        "bio/samtools/index"


rule samtools_cram:
    input:
        "star/{sample}/{maptype}/{sample}.bam",
        bam_index="star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
    output:
        protected("data_output/Mapping/{maptype}/{sample}.cram"),
    threads: min(config.get("max_threads", 2), 2)
    resources:
        mem_mb=get_6gb_per_attempt,
        time_min=get_2h_per_attempt,
        tmpdir="tmp",
    retries: 1
    params:
        extra=config["samtools"].get("gram_extra", ""),
    log:
        "logs/samtools/cram/{sample}.{maptype}.log",
    wrapper:
        "bio/samtools/view"
