rule bigr_copy:
    output:
        "reads/{status}/{sample}.{stream}.fq.gz"
    message:
        "Gathering {wildcards.status} {wildcards.sample} fastq files "
        "({wildcards.stream})"
    threads: 1
    resources:
        mem_mb=lambda wildcard, attempt: min(attempt * 1024, 2048),
        time_min=lambda wildcard, attempt: attempt * 45,
        tmpdir="tmp"
    params:
        input=lambda w, output: fastq_links[w.status][output[0]]
    log:
        "logs/bigr_copy/{status}/{sample}.{stream}.log"
    wrapper:
        "bio/BiGR/copy"