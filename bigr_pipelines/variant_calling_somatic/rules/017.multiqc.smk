rule multiqc:
    input:
        html=expand(
            "fastp/html/pe/{sample}_{status}.fastp.html",
            sample=design["Sample_id"],
            status=status_list,
        ),
        json=expand(
            "fastp/json/pe/{sample}_{status}.fastp.json",
            sample=design["Sample_id"],
            status=status_list,
        ),
        sambamba_metrics=expand(
            "sambamba/markdup/{sample}_{status}.bam",
            sample=design["Sample_id"],
            status=status_list,
        ),
        fastq_screen=expand(
            "fastq_screen/{sample}.{stream}.{status}.fastq_screen.{ext}",
            sample=design["Sample_id"],
            stream=streams,
            ext=["txt", "png"],
            status=status_list,
        ),
        picard=[
            multiext(
                f"picard/stats/{sample}_{status}.{cleaning}",
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
            for cleaning in cleaning_status
            for sample in sample_list
            for status in status_list
        ],
        csvstats=expand("snpeff/csvstats/{sample}.csv", sample=sample_list),
    output:
        protected("data_output/MultiQC/Somatic_Variant_Calling.html"),
    threads: 1
    resources:
        mem_mb=get_2gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    retries: 1
    log:
        "logs/multiqc.log",
    wrapper:
        "bio/multiqc"
