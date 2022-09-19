rule bigr_copy:
    output:
        "reads/{sample}.fq.gz"
    message:
        "Gathering {wildcards.sample} fastq file"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        input=lambda wildcards, output: fastq_links[output[0]]
    log:
        "logs/bigr_copy/{sample}.log"
    wrapper:
        str(worflow_source_dir / "bio" / "BiGR" / "copy")