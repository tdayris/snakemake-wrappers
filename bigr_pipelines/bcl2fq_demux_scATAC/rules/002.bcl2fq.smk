rule bcl2fastq:
    input:
        run_dir=config["run_dir"],
        sample_sheet=config["sample_sheet"],
        interop_dir=directory("InterOp")
    output:
        reports_dir=directory("unaligned/Reports/"),
        stats_json=directory("unaligned/Stats/")
        undetermined=expand(
            "unaligned/Undetermined_S0_{undetermined}_001.fastq.gz,
            undetermined=list_undetermined(config["dual"], config["umi"])
        )
    message:
        "Running bcl2fastq on {config['run_dir']}"
    threads: config.get("threads", 36)
    resources:
        mem_mb=lambda wildcards, attempt: min(attempt * 10240 + 20480, 51200),
        time_min=lambda wildcards, attempt: attempt * 60 * 12 * config.get("flowcell_size", 1)
    log:
        "logs/bcl2fastq/demux.log"
    params:
        extra=(
            "--fastq-compression-level 6 "
            "--mask-short-adapter-reads 1 "
            "--create-fastq-for-index-reads "
        ),
        out_dir=config["out_dir"],
        use_bases_mask=config.get("use_bases_mask", None)
        no_lane_splitting=config.get("no_lane_splitting", False),
        barcode_mismatches=config.get("barcode_mismatches", None)
    wrapper:
        "bio/bcl2fastq"