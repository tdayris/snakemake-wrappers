rule bigr_copy:
    output:
        temp("data_input/{sample_stream}.fq.gz")
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get45min_per_attempt,
        tmpdir="tmp",
    params:
        input=lambda wildcards, output: get_bigr(wildcards, output),
        extra="-v",
        extra_irods="-vK",
        extra_ln="-sfrv",
    log:
        "logs/020.bigr.copy/{sample_stream}.fq.gz"
    conda:
        str(workflow_source_dir / "envs" / "020.bigr.copy.yaml")
    script:
        str(workflow_source_dir / "scripts" / "020.bigr.copy.py")