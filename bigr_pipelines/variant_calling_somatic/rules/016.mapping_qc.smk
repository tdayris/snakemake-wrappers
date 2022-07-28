rule collect_multiple_metrics_raw:
    input:
        bam="bwa_mem2/sorted/{sample}_{status}.bam",
        bai="bwa_mem2/sorted/{sample}_{status}.bam.bai",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
    output:
        temp(
            multiext(
                "picard/stats/{sample}_{status}.raw",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".error_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
            )
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/picard/multiple_metrics/{sample}.raw.log",
    params:
        extra=config["picard"].get("collect_multiple_metrics", ""),
    wrapper:
        "bio/picard/collectmultiplemetrics"


rule collect_multiple_metrics_cleaned:
    input:
        bam="sambamba/markdup/{sample}_{status}.bam",
        bai="sambamba/markdup/{sample}_{status}.bam.bai",
        ref=config["reference"]["fasta"],
        ref_idx=config["reference"]["fasta_index"],
        ref_dict=config["reference"]["fasta_dict"],
    output:
        temp(
            multiext(
                "picard/stats/{sample}_{status}.cleaned",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
                ".bait_bias_detail_metrics",
                ".bait_bias_summary_metrics",
                ".error_summary_metrics",
                ".pre_adapter_detail_metrics",
                ".pre_adapter_summary_metrics",
            )
        ),
    threads: 1
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_1h_per_attempt,
        tmpdir="tmp",
    log:
        "logs/picard/multiple_metrics/{sample}.cleaned.log",
    params:
        extra=config["picard"].get("collect_multiple_metrics", ""),
    wrapper:
        "bio/picard/collectmultiplemetrics"


rule fastq_screen:
    input:
        "fastp/trimmed/{sample}_{status}.{stream}.fastq",
    output:
        txt=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.txt"),
        png=temp("fastq_screen/{sample}.{stream}.{status}.fastq_screen.png"),
    threads: config.get("threads", 20)
    resources:
        mem_mb=get_8gb_per_attempt,
        time_min=get_75min_per_attempt,
        tmpdir="tmp",
    params:
        fastq_screen_config=config["fastq_screen"],
        subset=100000,
        aligner="bowtie2",
    log:
        "logs/fastqc/{sample}.{stream}.{status}.log",
    wrapper:
        "bio/fastq_screen"
