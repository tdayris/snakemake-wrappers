rule sambamba_index_bam:
    input:
        "{tool}/{subcommand}/{sample}_{status}.bam",
    output:
        temp("{tool}/{subcommand}/{sample}_{status}.bam.bai"),
    threads: min(config.get("max_threads", 8), 8)
    resources:
        mem_mb=get_4gb_per_attempt,
        time_min=get_20min_per_attempt,
        tmpdir=tmp,
    retries: 1
    log:
        "sambamba/index/{tool}_{subcommand}/{sample}_{status}.log",
    params:
        extra=config["sambamba"].get("index_extra", ""),
    wrapper:
        "bio/sambamba/index"
