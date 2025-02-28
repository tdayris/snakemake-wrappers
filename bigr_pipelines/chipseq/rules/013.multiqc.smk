rule multiqc_report:
    input:
        picard=[
            multiext(
                f"picard/stats/{sample}",
                ".alignment_summary_metrics",
                ".insert_size_metrics",
                ".insert_size_histogram.pdf",
                ".quality_distribution_metrics",
                ".quality_distribution.pdf",
                ".gc_bias.detail_metrics",
                ".gc_bias.summary_metrics",
                ".gc_bias.pdf",
            )
            for sample in sample_list
        ],
        html=expand("fastp/html/{sample}.fastp.html", sample=sample_list),
        json=expand("fastp/json/{sample}.fastp.json", sample=sample_list),
        fingerprint_counts="deeptools/plot_fingerprint/raw_counts.tab",
        fingerprint_qc_metrics="deeptools/plot_fingerprint/qc_metrics.txt",
        correlations="deeptools/correlations/SpearmanCorr_readCounts.tab",
        pca="deeptools/pca/PCA.txt",
    output:
        "data_output/Report.html",
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir="tmp",
    log:
        "logs/multiqc.log",
    params:
        "",
    wrapper:
        "bio/multiqc"
