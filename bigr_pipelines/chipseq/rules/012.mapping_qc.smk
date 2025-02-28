rule collect_multiple_metrics:
    input:
        bam="bowtie2/sorted/{sample}.bam",
        bai="bowtie2/sorted/{sample}.bam.bai",
        ref=config["reference"]["genome"],
        ref_idx=config["reference"]["genome_index"],
        ref_dict=config["reference"]["genome_dict"],
    output:
        multiext(
            "picard/stats/{sample}",
            ".alignment_summary_metrics",
            ".insert_size_metrics",
            ".insert_size_histogram.pdf",
            ".quality_distribution_metrics",
            ".quality_distribution.pdf",
            ".gc_bias.detail_metrics",
            ".gc_bias.summary_metrics",
            ".gc_bias.pdf",
        ),
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/picard/multiple_metrics/{sample}.log",
    params:
        # optional parameters
        # REF_FLAT is required if RnaSeqMetrics are used
        extra=""
    wrapper:
        "bio/picard/collectmultiplemetrics"
