rule mixcr_plots_segment_usage:
    input:
        json="mixcr/post_analysis/{comparison}/individual_post_analysis.json.gz",
        metadata="mixcr/post_analysis/{comparison}/metadata.tsv"
    output:
        svg=protected("data_output/Mixcr/{comparison}/{segment}.svg")
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/mixcr/plots_segment_usage/{comparison}.{segment}.log"
    params:
        extra=config["mixcr"].get("plots_segment_usage", ""),
        subcommand=lambda wildcards: wildcards.segment
    wrapper:
        "bio/mixcr/exportPlots"


rule mixcr_plot_overlap:
    input:
        json="mixcr/post_analysis/{comparison}/overlap_post_analysis.json.gz",
        metadata="mixcr/post_analysis/{comparison}/metadata.tsv"
    output:
        svg=protected("data_output/Mixcr/{comparison}/Overlap.svg")
    threads: config.get("max_threads", 20)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_3h_per_attempt,
        tmpdir="tmp"
    log:
        "logs/mixcr/plots_segment_usage/{comparison}.log"
    params:
        extra=config["mixcr"].get("plots_overlap", ""),
        subcommand="overlap"
    wrapper:
        "bio/mixcr/exportPlots"