rule collect_multiple_metrics:
    input:
        bam="star/{sample}/{maptype}/{sample}.bam",
        bai="star/{sample}/{maptype}/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
        refflat=config["reference"]["refflat"],
    output:
        temp(
            multiext(
                "picard/stats/{sample}.{maptype}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
                ".rna_metrics",
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".error_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
            )
        ),
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/picard/multiple_metrics/{sample}.{maptype}.log",
    params:
        extra=lambda wildcards, input: f"REF_FLAT {input.refflat} {config['picard'].get('extra', '')}",
    wrapper:
        "bio/picard/collectmultiplemetrics"


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
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp",
    retries: 2
    params:
        extra=config["samtools"].get("gram_extra", ""),
    log:
        "logs/samtools/cram/{sample}.{maptype}.log",
    wrapper:
        "bio/samtools/view"
