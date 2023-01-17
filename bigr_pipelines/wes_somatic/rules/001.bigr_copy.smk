rule bigr_copy:
    output:
        temp("data_input/{status}/{sample}.{stream}.fq.gz"),
    threads: 1
    resources:
        mem_mb=get_1gb_per_attempt,
        time_min=get_45min_per_attempt,
        tmpdir=tmp,
    retries: 1
    params:
        input=lambda w, output: fastq_links[w.status][output[0]],
    log:
        "logs/bigr_copy/{status}/{sample}.{stream}.log",
    wrapper:
        "bio/BiGR/copy"
