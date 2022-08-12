rule material_methods:
    output:
        protected("data_output/Mat&Met.txt")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_15min_per_attempt,
        tmpdir="tmp"
    log:
        "logs/mat&met.log"
    params:
        all_parameters=config.copy()
    script:
        str(workflow_source_dir / "scripts" / "mat_n_met.py")